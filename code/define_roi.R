# Define ROI and land line layers
library(tidyverse)
library(sf)
library(terra)
library(rworldxtra)
data("countriesHigh")
library(elevatr)

# Load AOI layers
world = countriesHigh %>% st_as_sf() %>% st_transform(32616)
land = world %>% filter(REGION %in% c('North America', 'South America and the Caribbean'))
sf_use_s2(TRUE)

# extent of interest
ylim = c(5,50)
xlim = c(-120,-60)
roi = tibble(lat = ylim[c(1,2,2,1,1)],
             lon = xlim[c(1,1,2,2,1)]) %>% 
  st_as_sf(coords = c('lon', 'lat'), crs=4326) %>%
  st_bbox() %>%
  st_as_sfc() %>%
  st_transform(32616) %>%
  st_bbox() %>%
  st_as_sfc()

land_buff_200km = land %>% st_buffer(dist=1) %>% st_union() %>% st_buffer(dist=200000)

ggplot() + 
  geom_sf(data=land_buff_200km, linetype='dashed', size=0.1, fill = 'black', alpha=0.1) +
  geom_sf(data=land, mapping = aes(fill=ne_10m_adm), size=0.1) +
  theme_bw() + theme(legend.position='none') +
  geom_sf(data=roi, linetype='dashed', fill=NA) +
  lims(x=st_bbox(roi)[c(1,3)], y = st_bbox(roi)[c(2,4)]) +
  scale_fill_viridis_d(option='D') +
  theme(panel.background = element_rect(fill = 'lightblue'))
  
ggsave('figs/roi.png', width=5, height=4, dpi=600)

land_buff_200km = land_buff_200km %>% st_crop(roi)  
land = land  %>% st_crop(roi) %>% st_union()

## OUTPUT land and ROI shapefiels
write_sf(land, 'input/land.shp')
write_sf(roi, 'input/roi.shp')
write_sf(land_buff_200km, 'input/land_buff200km.shp')
write_sf(world %>% filter(REGION %in% c('North America', 'South America and the Caribbean')), 'input/land_lines.shp')
