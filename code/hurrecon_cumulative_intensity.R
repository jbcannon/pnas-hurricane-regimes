setwd('/home/jcannon/hurrecon-rev')

# Script to combine all hurrecon outputs into rasters of probablities of windspeed
# at different intensities SS 0- SS5

wind_breaks_mph = c(39, 74, 95, 110, 129, 156, 2000)
res_m = 1000
path_to_hurrecon_Vs_outputs = 'data/'
n_workers = 3 #Processes RAM intensive. On 126GB mem, 4 workers unstable
wind_breaks = wind_breaks_mph * 0.44704 #convert mph to m/s

# Load libraries
library(parallel)
library(foreach)
library(doParallel)
library(terra)
library(sf)
library(hurrecon)
library(ggplot2)

#Load ROI info
roi = read_sf('input/roi.shp')
land = read_sf('input/land.shp')
land_lines = read_sf('input/land_lines.shp')
land_buff_200km = read_sf('input/land_buff200km.shp')

#Function to divide parallel tasks among workers
divide_tasks_among_workers = function(task_iterator, n_workers, shuffle=FALSE) {
  n_tasks = length(task_iterator)
  factor = floor(n_tasks/n_workers)
  remain = n_tasks %% (factor * n_workers)
  if(shuffle) task_iterator = sample(task_iterator)
  assignments = list()
  for(i in 1:n_workers) assignments[[i]] = task_iterator[(1:factor)+(i-1)*factor]
  if (remain > 0) for(j in 1:remain) assignments[[j]] = c(assignments[[j]], task_iterator[(n_workers*factor + j)])
  return(assignments)
}

divide_into_chunks = function(tasks, chunks, shuffle=TRUE) {
  output = list()
  if(shuffle) tasks = sample(tasks)
  i=1
  while(length(tasks) > 0) {
    if(length(tasks) >= chunks) {
      output[[i]] =  tasks[1:chunks]
      tasks = tasks[-(1:chunks)]
    } else {
      output[[i]] = tasks
      tasks = NULL
    }
   i=i+1
  }
  names(output) = 1:length(output)
  return(output)
}

# Get list of all files to combine
trks = list.files(path_to_hurrecon_Vs_outputs, pattern = 'Vs.tif', full.names = TRUE)

#Check all these tracks to make sure they only have 1 layer
# an earlier bug was creating multi-layer rasters that slowed and crashed things
#for(i in seq_along(trks)) if(terra::nlyr(rast(trks[i]))>1) stop('yep')

# Break apart EP and AL hurricanes to divide by appropriate time seriess
al_trks = trks[substr(trks,7,8)=='AL']
ep_trks = trks[substr(trks,7,8)=='EP']
al_assignments = divide_into_chunks(al_trks, 10)
ep_assignments = divide_into_chunks(ep_trks, 10)
assignments = c(al_assignments, ep_assignments)
names(assignments) = seq_along(assignments)

unlink('cluster-log-cumulative-intensity.txt')
registerDoParallel(makeCluster(n_workers, outfile = 'cluster-log-cumulative-intensity.txt'))

#Create list of frequency rasters
SS_out =  foreach(
 #j = seq_along(assignments),
  j= 1:length(assignments),
 .packages = c('hurrecon', 'terra', 'sf'),
 .combine = c,
 .init = c()
) %dopar% {
  basin_name = substr(assignments[[j]][1],7,8)
  chunk_name = paste0(basin_name, '_', sprintf('%03d', j))
  counter = 1
  cat(paste0(chunk_name,' started\n'))

# Create raster template within each worker (large rast is non-exportable)
  e = st_bbox(roi)
  e[c(1,3)] = floor((e[c(1,3)]/res_m))*res_m
  e[c(2,4)] = ceiling((e[c(2,4)]/res_m))*res_m
  x = rast(xmin=e[1], xmax=e[3], ymin=e[2], ymax=e[4],
           crs = st_crs(roi)$wkt)
  res(x) = rep(res_m,2)
  x[] = 0
  
  # Run parallel combination script
  i=assignments[[j]]
  SS_freq = c(x,x,x,x,x,x)
  for(t in i){
    tmp_surf = rast(t)
    for(ss in 5:0) if(sum((tmp_surf>=wind_breaks[ss+1])[])) {ss_list = 0:ss; break}
    tmp_surf = terra::crop(tmp_surf, x)
    tmp_surf = terra::extend(tmp_surf, x)
    tmp_surf[is.na(tmp_surf)] = 0
    for(ss in ss_list) {
      tmp_ss = tmp_surf >= wind_breaks[ss+1] &  tmp_surf < wind_breaks[ss+2]
      SS_freq[[ss+1]] = SS_freq[[ss+1]] + tmp_ss
      cat(chunk_name, '-', counter, ':\t', basename(t), 'SS:', ss, 'completed\n')
    }
    counter=counter+1
  }
  # Write to temp file
  fn = paste0(getwd(), '/tmp/assgn_', chunk_name, '.tif')
  writeRaster(SS_freq, fn)
  gc()
  return()
}

# Calculate number of instances of SS for each basin
for(b in c('EP', 'AL')) {
  SS_out = list.files('tmp/', pattern=paste0('assgn_',b), full.names=TRUE)
  x = lapply(SS_out, function(x) rast(x))
  y = x[[1]]
  for(i in 2:length(SS_out)) {y = y + x[[i]]; cat(i,'\n')}
  SS_final = y
  names(SS_final) = paste0('ss_', 0:5)
  writeRaster(SS_final, paste0('output/SS_freq_tally-', b,'.tif'), overwrite=TRUE)
}

## Pull both basins in and get risk..
SS_ep = rast('output/SS_freq_tally-EP.tif')
SS_al = rast('output/SS_freq_tally-AL.tif')

df = tibble(fn = trks,
            id = gsub('data//|_Vs.tif', '',fn),
            yr = substr(id, 5,8),
            basin = substr(id,1,2))
df %>% group_by(basin) %>% summarize(yr_range = range(yr))
yr_rng_al = 2021-1851+1
yr_rng_ep = 2022-1949+1

risk_ep = (SS_ep/yr_rng_ep)
risk_al = (SS_al/yr_rng_al)

# Calculate 1-yr 10-yr risk and 30-yr risk from hurricane tallys
risk_ep_1yr = aggregate(risk_ep,5,fun=mean) %>%
  focal(w=matrix(1,5,5,5), fun=mean)
risk_al_1yr = aggregate(risk_al,5,fun=mean) %>%
  focal(w=matrix(1,5,5,5), fun=mean)

risk_1yr = risk_ep_1yr + risk_al_1yr
risk_10yr = rast(lapply(1:nlyr(risk_1yr), function(i) 1-(1-risk_1yr[[i]])^10))
risk_30yr = rast(lapply(1:nlyr(risk_1yr), function(i) 1-(1-risk_1yr[[i]])^30))

#Write risk layers to disk
writeRaster(risk_1yr, 'output/SS_risk_01yr.tif', overwrite=TRUE)
writeRaster(risk_10yr, 'output/SS_risk_10yr.tif', overwrite=TRUE)
writeRaster(risk_30yr, 'output/SS_risk_30yr.tif', overwrite=TRUE)

####  Plot outputs

# Generate plot of 1, 10, and 30 yr risk
its = c('01', '10', '30') 
for(it in its) {
  png(paste0('figs/SS_',it,'yr_risk_panel.png'), res=300, width=12, height=7, units='in')
  risk = rast(paste0('output/SS_risk_',it,'yr.tif'))
  par(mfrow=2:3, oma=c(0,0,2,1.25))
  for(i in 1:6){
    rng = range(risk[[i]][], na.rm=TRUE)
    ras_breaks_tmp = pretty(c(0, rng[2]), n=5)
    #plot(roi$geometry, asp=1, axes=FALSE)
    terra::plot(risk[[i]], breaks = ras_breaks_tmp, col=topo.colors(length(ras_breaks_tmp)-1), axes=FALSE, extent=roi, mar=c(0,0,0,0))
    plot(st_crop(land_lines,roi)$geometry,add=1)
    mtext(paste0('Saffir-Simpson\nCategory ', i-1), 3,-3.5,cex=0.7,adj=0.15, col='white')
  }
  mtext(paste0(as.numeric(it),' year hurricane risk'), 3, outer=TRUE)
  dev.off()
}
