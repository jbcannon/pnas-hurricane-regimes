#setwd('/home/jeffery.cannon@jonesctr.org/hurrecon-rev')

#download latest version of hurrecon
devtools::install_github('jbcannon/hurrecon')
library(hurrecon)
data('geographic')
data('radius_models')
library(tidyverse)
library(sf)
library(terra)
library(doParallel)
library(parallel)
library(foreach)

# Load hurrdat data
hurdat = 'input/AL-NECP-combined.csv'

# Output folder
output_folder = 'data/'

# Number of cores
useCores = 15

land_buff_200km = read_sf('input/land_buff200km.shp')
roi = read_sf('input/roi.shp')

# Select which tracks to work on
if(!file.exists(hurdat)) hurrecon::fetch_best_tracks_data(hurdat)

SS_summary = read_csv(hurdat) %>% filter(status=='HU') %>%
  mutate(ss = cut(max_speed, c(0,64,82,96,112,137,500), labels = 0:5) %>% as.numeric()) %>%
  group_by(track_id) %>% summarise(SS = max(ss))
track_list = SS_summary$track_id
done = output_folder %>% list.files('_Vs.tif|.csv') %>% str_replace('_Vs.tif|.csv|ERROR_', '') %>% str_replace('_Vs.tif|.csv','')
todo_list = track_list[!track_list %in% done]

# order by year to give idea of progress
todo_list = todo_list[substr(todo_list, 5,8) %>% as.numeric() %>% order()]

# See if there are any old error files (remove comment to delete them)
output_folder %>% list.files(pattern='ERROR', full.names=TRUE)  #%>% unlink()

# Continue to run in parallel or skip below L111 to run in series.

# Parallel options
useCores %>% makeCluster(outfile='cluster-log.txt') %>% registerDoParallel

packages = c('hurrecon', 'terra', 'tidyverse')

print(todo_list)
#--> Loop in parallel
par_out = foreach(
  i = todo_list,
  .packages = packages,
  .errorhandling = 'remove') %dopar% {
    data("radius_models")
    
    trk = load_hurdat_track(hurdat, i)
    
    # Check to see if its already been done
    done_trks = output_folder %>% list.files('_Vs.tif|.csv') %>% str_replace('_Vs.tif|.csv|ERROR_', '') %>% str_replace('_Vs.tif|.csv','')
    done = i%in%done_trks
    
    # Check to see if it touches land in the the area of interest
    touch_land = any(as.numeric(sf::st_intersects(trk, land_buff_200km, sparse = TRUE)) == 1)
    
    #check to see if it is in the roi
    trk = trk[land_buff_200km,]
    if(nrow(trk)>=2) trk = trk[roi,]
    in_roi = nrow(trk)>=2
    
    #Check to see if it becomes a hurricane
    hurr = any(trk$status=='HU')
    if(!done & hurr & touch_land & in_roi) {
      x = try(hurrecon_run(trk, max_rad_km = 300, res_m = 1000, max_interp_dist_km = 2, aoi=roi, land =land))
      if('try-error' %in% class(x)) {
        err_tab = data.frame(track_id=i,note='error')
        readr::write_csv(err_tab, file=paste0(output_folder, '/ERROR_', i,'.csv'))
      } else {
        writeRaster(x, paste0(output_folder, '/', i, '_Vs.tif'), overwrite=TRUE)
      }
    } else {
      if(!done) {
        err_tab = data.frame(track_id=i,note='invalid')
        readr::write_csv(err_tab, file=paste0(output_folder, '/', i,'.csv'))
      }
    }
  }

## Summary of how many done
n = 'data/' %>% list.files('.tif') %>% length

rg = 'data/' %>% list.files('.tif') %>% substr(5,8) %>% 
  as.numeric() %>% range

cat(n, 'events identified between', rg[1], 'and', rg[2],
    'corresponding to ', round(n/diff(rg),2), 'per year')


### RUN in series


data("radius_models")


# Find tracks that are not yet complete
done = output_folder %>% list.files('_Vs.tif|.csv') %>% str_replace('_Vs.tif|.csv|ERROR_', '') %>% str_replace('_Vs.tif|.csv','')
todo_list = track_list[!track_list %in% done] %>% rev

# OR find tracks that ran errors
remain = output_folder %>% list.files('ERROR') %>% str_replace('_Vs.tif|.csv|ERROR_', '') %>% str_replace('_Vs.tif|.csv','')
todo_list = track_list[track_list %in% remain] %>% rev

# Go through each individual (or add a loop if you want, but this is for troubleshooting)

i = todo_list[3]

trk = load_hurdat_track(hurdat, i)

# Check to see if it touches land in the the area of interest
touch_land = any(as.numeric(sf::st_intersects(trk, land_buff_200km, sparse = TRUE)) == 1)

#check to see if it is in the roi
trk = trk[land_buff_200km,]
if(nrow(trk)>=2) trk = trk[roi,]
in_roi = nrow(trk)>=2

#Check to see if it becomes a hurricane
hurr = any(trk$status=='HU')
if(hurr & touch_land & in_roi) {
  x = try(hurrecon_run(trk, max_rad_km = 300, res_m = 1000, max_interp_dist_km = 2, aoi=roi, land=land))
  if('try-error' %in% class(x)) {
    err_tab = data.frame(track_id=i,note='error')
    readr::write_csv(err_tab, file=paste0(output_folder, '/ERROR_', i,'.csv'))
  } else {
    writeRaster(x, paste0(output_folder, '/', i, '_Vs.tif'), overwrite=TRUE)
  }
} else {
  if(!done) {
    err_tab = data.frame(track_id=i,note='invalid')
    readr::write_csv(err_tab, file=paste0(output_folder, '/', i,'.csv'))
  }
}

## How many were modeled in each basin
x = list.files(output_folder, '.Vs.tif')
df = tibble(fn=x, basin=substr(fn,1,2))
table(df$basin)
