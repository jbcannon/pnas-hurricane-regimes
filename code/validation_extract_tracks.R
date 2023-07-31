#------------------------------------------------------------------------------>
# IDENTIFY HURRICANES TO VALIDATE WITH HRR
#------------------------------------------------------------------------------>
# 
# Jeffery B. Cannon
# 23 June 2022
# The Jones Center at Ichauway
# Landscape Ecology Lab
# jeffery.cannon@jonesctr.org
#
# This script selects hurricane tracks from Hurrecon study that
#  - occured on or after 2014
#  - stratified among category I - V storms
#  - spent ~ 3000 km on land
#
#------------------------------------------------------------------------------>
# LOAD LIBRARIES AND DATA
#------------------------------------------------------------------------------>

library(ggplot2)
library(tidyverse)
library(dplyr)
#install_github('jbcannon/hurrecon')
library(hurrecon)
data(geographic)
hurdat = 'input/AL-NECP-combined.csv'

#------------------------------------------------------------------------------>
# PROCESSING
#------------------------------------------------------------------------------>

# Subset hurricane tracks to only post 2014 and in Atlantic basin
if(!file.exists(hurdat)) hurrecon::fetch_best_tracks_data(hurdat)
tracks = read_csv(hurdat) %>% group_by(track_id) %>%
  summarise(HU = any(status=='HU'),
            year = str_extract(date, "^.{4}")[1]) %>%
  filter(year >= 2014 & HU == TRUE) %>%
  filter(substr(track_id,1,2) == 'AL') %>%
  pull(track_id)

# Convert track points into line objects
track_lines = read_csv(hurdat) %>% filter(track_id %in% tracks) %>%
  st_as_sf(coords = c('lon', 'lat'))
out = list()
for(i in unique(track_lines$track_id)) {
  tmp = filter(track_lines, track_id == i) %>% 
    summarise(do_union=FALSE) %>% st_cast('LINESTRING')
  st_crs(tmp) = 4326
  tmp$track_id = i
  out[[length(out)+1]] = tmp
}
track_lines = do.call(rbind, out)

# Summarize track-level info (start/end time)
track_info_full = read_csv(hurdat) %>% filter(track_id %in% tracks) %>%
  mutate(datetime = paste0(date,time) %>% as.numeric()) %>%
  group_by(track_id) %>% summarise(vmax=max(max_speed),
                                   full_s = min(datetime),
                                   full_e = max(datetime))
# Which lines "touch" land
temp_land = st_transform(land, 32616)
track_lines = st_transform(track_lines, 32616)
track_lines = track_lines[st_intersects(track_lines, temp_land, sparse = FALSE),]
tracks = track_lines$track_id

# What storms come withing 250 km of land?
land_buff = st_buffer(land, dist=250000) %>% st_transform(st_crs(track_lines))

# Get start/end times of approaching land
track_pts = read_csv(hurdat) %>% filter(track_id %in% tracks) %>%
  st_as_sf(coords = c('lon', 'lat'), crs=4326) %>%
  st_transform(st_crs(land_buff)) %>% 
  st_intersection(land_buff)

track_info = track_pts %>%
  group_by(track_id) %>%
  mutate(datetime = paste0(date,time) %>% as.numeric()) %>%
  summarise(vmax = max(max_speed),
            int_s = min(datetime),
            int_e = max(datetime)) %>% 
  st_drop_geometry()

out = list()
for(i in unique(track_pts$track_id)) {
  tmp = filter(track_pts, track_id == i) %>% 
    summarise(do_union=FALSE) %>% st_cast('LINESTRING')
  #st_crs(tmp) = 4326
  tmp$track_id = i
  out[[length(out)+1]] = tmp
}
track_pts = do.call(rbind, out)

track_ints = left_join(track_pts, track_info) %>%
  left_join(track_info_full)
plot(track_ints)

#------------------------------------------------------------------------------>
# OUTPUT DATA
#------------------------------------------------------------------------------>

st_write(track_ints, 'validation/track_segments.shp', append = FALSE)

track_ints = st_read('validation/track_segments.shp')

files = grep(pattern = track_ints$track_id %>% paste0(collapse='|'), 
     x= list.files('data/', pattern='.tif', full.names=TRUE), value = TRUE)
zip('validation/hurrecon_outputs.zip', files)

#------------------------------------------------------------------------------>
# PLOT
#------------------------------------------------------------------------------>

xlim = st_bbox(land)[c(1,3)]
xlim = c(-1e6,3e6)
ylim = st_bbox(land)[c(2,4)]
ylim = c(11e5,5e6)

plot_out = ggplot(land) + geom_sf(fill=grey(0.8)) +
  geom_sf(data = land_lines, fill=NA) +
  geom_sf(data = track_ints, lineend='round', alpha=0.2, aes(size=vmax)) +
  geom_sf(data=land_buff, fill=NA, linetype='dashed') +
  coord_sf(x=xlim, y=ylim)

png('validation/included_tracks.png', width=6, height=6, res=600, units='in')
plot(plot_out)
dev.off()