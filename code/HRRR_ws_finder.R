library(sf)
library(terra)
library(tidyverse)
library(hurrecon)

hrrr_to_raster = function(hrr_fn, lat_ras, lon_ras, out_file=NULL)
  {
  if(!'lon' %in% ls()) {
  lon <<- ((read_table('HRRR_outputs/longitude.txt', col_names=FALSE)  %>%
             as.matrix() %>% rast) -360) %>% flip}
  if(!'lat' %in% ls()) {
  lat <<- read_table('HRRR_outputs/latitude.txt', col_names=FALSE) %>%
      as.matrix() %>% rast %>% flip}
  
  ws = read.table(hrr_fn) %>% as.matrix %>% rast
    df = tibble(ws=as.numeric(ws[]),
                lon=as.numeric(lon[]),
                lat=as.numeric(lat[]))
    sp = st_as_sf(df, coords=c('lon', 'lat'))
    st_crs(sp) = st_crs(4326)
    r = rast(ncols = ncol(ws), nrows=nrow(ws))
    ext(r) = ext(sp)
    x = terra::rasterize(vect(sp), r, fun=mean, field='ws') %>%
      focal(w=3, fun=mean, na.policy='only', na.rm=T)
    if(is.null(out_file)) return(x)
    writeRaster(x, out_file)
    return(out_file)
}

# Run this for all the HRR outputs
hrr_list = list.files('HRRR_outputs/', pattern = '.npy.txt', full.names = TRUE)

for(i in hrr_list) {
  cat(i)
  cat('\nfile', which(i==hrr_list), 'of', length(hrr_list),'\n')
  out_file = gsub('.npy.txt', '.tif', i)
  if(file.exists(out_file)) next
  x = hrrr_to_raster(i,
                 lat_ras=lat,
                 lon_ras=lon,
                 out_file = gsub('.npy.txt', '.tif', i))
  }


### Load test tracks and compare them

tracks = read_sf('track_segments.shp') %>% st_transform(st_crs(4326))
data('geographic')

# Make polygon with hrrr extent
x = rast('HRRR_outputs/WS_2014-08-01.tif')
x[!is.na(x)] = 0
x = as.polygons(x) %>% st_as_sf()
hrr_ext = x; rm(x)
# get polygon of land cover
land = land %>% st_transform(st_crs(4326))

#Function to stratify sample by distance
#runif(1,min=0,max=1e6)
set.seed(179613)

spatial_strat_sample = function(f, breaks, n) {
  if(n %% length(breaks) != 0) stop('n must be a multiple of the number of breaks')
  orig_proj = st_crs(f)
  f = st_transform(f, st_crs(32616))
  f = lapply(breaks, function(i) st_buffer(f, dist=i))
  samps = n / length(breaks)
  pts = list()
  for(i in 1:length(f)) {
    shape = f[[i]]
    if(i > 1) shape = st_difference(shape, f[i-1][[1]])
    pts[[length(pts) + 1]] = st_sample(shape, samps)
  }; pts = do.call(c, pts)
  pts = st_transform(pts, orig_proj)
  pts = st_as_sf(pts)
  return(pts)
}

# Loop through all tracks and pull together
# HRR and VS predictions

full_df = list()
for(i in 1:nrow(tracks)) {
  # Load next track
  t = tracks[i,] 
  trk_name = t$track_id
  trk_date = t$int_s %>% substr(1,8) %>% as.Date('%Y%m%d')
  cat('working on ', trk_name, '\ntrack', i, 'of', nrow(tracks), '\n')
  
  # check if vs and hrr are known
  hrr_fn = paste0('HRRR_outputs/WS_', trk_date, '.tif')
  if(!file.exists(hrr_fn)) next
  vs_fn = paste0('hurrecon_outputs/', trk_name, '_Vs.tif')
  if(!file.exists(vs_fn)) next
  
  # Load Vs and HRR rasters
  vs = rast(vs_fn)
  hrrr = rast(hrr_fn)
  vs = project(vs, hrrr)
  names(vs) = 'vs'
  names(hrrr) = 'hrrr'
  
  # Crop layers to smaller extent
  t = st_intersection(t, hrr_ext)
  if(nrow(t) == 0) next
  t_buff = t %>% st_transform(st_crs(32616)) %>%
    st_buffer(dist=200000) %>% st_transform(st_crs(4326))
  hrrr = crop(hrrr, t_buff)
  vs = crop(vs, t_buff)
  
  # Select random points to compare among hrr and vs
  samp_pts = spatial_strat_sample(t,
                                  breaks=c(50,100,150,200)*1000,
                                  n=100)
  # Sample points and compare
  df = terra::extract(hrrr, vect(samp_pts)) %>%
    cbind(terra::extract(vs, vect(samp_pts))) %>%
    select(hrrr,vs) %>%
    mutate(id = trk_name)
  full_df[[length(full_df) + 1]] = df
  
  # Make a map
  png(paste0('validation_figures/', trk_name, '.png'), width=6, height=4, units='in', res=600)
  lims = range(c(hrrr[], vs[]), na.rm=TRUE) %>% pretty(10)
  par(mfrow=1:2, oma=c(0.5,1,2,0))
  plot(vs, mar=c(1.5,1,0,2), breaks=lims, legend=FALSE);
  plot(t$geometry, add=TRUE)
  plot(t_buff$geometry,add=TRUE, lty=2)
  points(t$geometry %>% st_coordinates, pch=16, cex=0.5)
  plot(samp_pts, add=TRUE, pch=3, cex=0.2, col='red')
  mtext('a. HURRECON', side=3,adj=0.95,line=-1.25, cex=0.8)
  plot(hrrr,mar=c(1.5,0,0,4), breaks=lims)
  plot(t$geometry, add=TRUE)
  points(t$geometry %>% st_coordinates, pch=16, cex=0.5)
  mtext(trk_name, side=3, line=0.5, outer=TRUE)
  mtext('b. HRRR', side=3,adj=0.95,line=-1.25, cex=0.8)
  mtext('Latitude', side=2,line=0,outer=TRUE, cex=0.8)
  mtext('Longitude', side=1,line=-0.5,outer=TRUE, cex=0.8)
  plot(t_buff$geometry,add=TRUE, lty=2)
  plot(samp_pts, add=TRUE, pch=3, cex=0.2, col='red')
  dev.off()
}

full_df = do.call(rbind, full_df)
#write_csv(full_df, 'validation_output.csv')
full_df = read_csv('validation_output.csv')
full_df = subset(full_df, !(is.na(full_df$hrrr) | is.na(full_df$vs)))

hurrecon_validation_mod = lm(vs~hrrr, data=full_df)
h = hurrecon_validation_mod
#calculate RMSE
RSS = c(crossprod(h$residuals))
MSE = RSS / length(h$residuals)
RMSE = sqrt(MSE)

r2 = summary(hurrecon_validation_mod)$r.squared %>% round(3)
r2 = bquote(R^2==.(r2))
pval = (hurrecon_validation_mod %>% summary)[[4]][2,4]
pval = ifelse(pval < 0.0001, "P < 0.0001", paste0('P = ', round(pval,3)))
rmse = paste0('RMSE = ', round(RMSE,2))

validation_fig = full_df %>% ggplot(aes(x=hrrr, y=vs)) +
  geom_point(alpha=0.4) +
  geom_smooth(method='lm') +
  geom_abline(slope=1,intercept=0, linetype='dashed') +
  coord_fixed() + 
  theme_bw() +
  labs(x = bquote(v[s~HRRR]~~'('*m~s^-1*')'),
       y = bquote(v[s~HURRECON]~~'('*m~s^-1*')')) +
  annotate('text', x=10,y=60, label=pval, size=3) +
  annotate('text', x=10,y=55, label=r2, size=3) +
  annotate('text', x=10,y=65, label=rmse, size=3)

print(validation_fig)

png('validation-fig.png', width=3, height=4, units='in', res=600)
print(validation_fig)
dev.off()

sink('validation-results.txt')
cat(rep('#', 20))
cat('\nHURRECON vs. HRR validation\n')
cat(rep('#', 20))
cat('\n\n')
cat('Number of tracks:', length(unique(full_df$id)))
cat('\n\n~ 100 points per track stratified by distance from track\n')
print(summary(hurrecon_validation_mod))
cat('\n\n')
cat('Tracks included\n\n')
incl = tracks %>% 
  filter(track_id %in% full_df$id) %>%
  select(track_id, vmax, int_s) %>%
  st_drop_geometry()
print(incl, n = nrow(incl))
sink()
