### Cluster analysis for 30 yr hurricane risk
library(terra)
library(tidyterra)
library(hurrecon)
library(RColorBrewer)
library(tidyverse)
library(lubridate)
library(ggplot2)
library(ggpubr)
library(ggridges)
library(ggnewscale)
library(grid)
library(rworldxtra)
data("countriesHigh")

# Load AOI layers
world = countriesHigh %>% st_as_sf() %>% st_transform(32616)
world = world %>% filter(REGION %in% c('North America', 'South America and the Caribbean'))
sf_use_s2(FALSE)

data('geographic')

# Load risk maps
riskmap_30 = rast('output/SS_risk_30yr.tif')
riskmap_10 = rast('output/SS_risk_10yr.tif')
riskmap_01 = rast('output/SS_risk_01yr.tif')

# Find areas with NO hurricane winds recorded
mask = lapply(1:nlyr(riskmap_30), function(ss) {
  tmp = riskmap_30[[ss]]
  tmp[tmp == 0] = NA
  return(tmp)}) %>% rast
mask = lapply(1:nlyr(mask), function(ss) return(1-is.na(mask[[ss]]))) %>% rast
mask = sum(mask)>0
mask[mask==0] = NA

x = fillHoles(as.polygons(mask)) 
mask = rasterize(x, mask)
x = mask(mask, land)
plot(x)
mask = x

# Delete areas where no hurricane recorded--> set to NA
riskmap_30_masked = lapply(1:nlyr(riskmap_30), function(ss) mask(riskmap_30[[ss]], mask))
riskmap_30 = rast(riskmap_30_masked)
z_score = function(x) (x - mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE)

# generate matrix with z-scores of 30yr hurricane risk for each category
x = riskmap_30
for(ss in 1:nlyr(x)) x[[ss]][] = z_score(x[[ss]][])
z_riskmap_30 = x

len = length(z_riskmap_30[[1]][])
mat = matrix(ncol= nlyr(z_riskmap_30), nrow = len)
for(ss in 1:nlyr(z_riskmap_30)) mat[,ss] = z_riskmap_30[[ss]][]
na_vals = is.na(mat[,1])
mat_slim = mat[!na_vals,]

# k-means cluster analysis of hurricane risk probabilities

# Scree plot to determine number of dimensions to use
wss <- (nrow(mat_slim)-1)*sum(apply(mat_slim,2,var))
for (i in 2:20) wss[i] <- sum(kmeans(mat_slim, centers=i, iter.max=100)$withinss)
wss = tibble(wss=wss, n = 1:20)

ggplot(data=wss, aes(x=n, y=wss)) + geom_point() + geom_line() +
  theme_classic() + labs(x='clusters', y = bquote(SS[within]))
ggsave('figs/cluster-analysis-scree-plot.png', width=6, height=6, dpi=600)

n_means = 4
set.seed(42665)
km = kmeans(mat_slim, n_means, iter.max=100)

#generate raster with clusters displayed
clus_ras = riskmap_30[[1]]
clus_ras[] = n_means+10
clus_ras[na_vals] = NA
clus_ras[!na_vals] = km$cluster

#Clusters don't necessarily come out ranked, so rank order them based on 30 yr risk (SS IV)
out = list()
out = do.call(rbind, lapply(1:nlyr(riskmap_30), function(ss) {
  tmp_out = zonal(riskmap_30[[ss]], clus_ras, fun =mean)
  tmp_out$ss = ss - 1
  colnames(tmp_out) = c('cluster', 'risk', 'ss')
  return(tmp_out)
}))
risk_summary = out

ggplot(data=out, aes(y=risk, x=cluster, group=ss, fill=ss)) +
  geom_bar(position='dodge', stat='identity')

# Make up some kind of cumulative expected wind speed (^2) to make it similar
# to energy (ie, like ACE index).
wind_breaks = tibble(ss=0:5, ws = c(17,33,42,49,58,70), ws2 = ws^2)
out = left_join(out, wind_breaks)
out$exp = out$risk * out$ws2
out = out %>% group_by(cluster) %>% summarize(cumulative_ex_ws = sum(exp))

# Use this cumulative index to re-order/rank by severity
sev_order = rank(out$cumulative_ex_ws)
out$cluster = as.factor(out$cluster)
out$cluster = factor(out$cluster, levels(out$cluster)[order(sev_order)])
levels(out$cluster) = 1:n_means

risk_summary$cluster = as.factor(risk_summary$cluster)
levels(risk_summary$cluster) = levels(risk_summary$cluster)[sev_order] 
risk_summary$cluster=as.factor(as.numeric(as.character(risk_summary$cluster)))
ggplot(data=risk_summary, aes(y=risk, x=cluster, group=ss, fill=ss)) +
  geom_bar(position='dodge', stat='identity')

# Now reassign values based on new names
tmp = clus_ras + 10
for(i in 1:n_means) {
  val = i + 10
  new_val = out$cluster[i] #get new renumbered value
  tmp[tmp==val] = new_val
}
clus_ras = tmp

#and create the vector version
x = as.polygons(clus_ras)
x$cluster = 1:n_means
x = disagg(x)
x = st_as_sf(x)
x$area_m2 = as.numeric(st_area(x))
x$area_km2 = x$area_m2/(1000*1000)
x$area_m2 = NULL
x[,1] = NULL
x = filter(x, area_km2 < 2.06e5 | area_km2 >2.07e5) #eliminate weird great lakes thing
clus_ras_vect = x

# Name these things
par(mfrow=c(2,2), mar=c(0,0,2,0))
titles = c('Continental', 'Inland', 'Coastal', 'Fringe')
for(i in 1:n_means) {
  col = 'cyan'; border=NA
  #if(i>3) {col='cyan'; border = NA}
  plot(filter(clus_ras_vect,cluster==i)$geom, col=col, border=border,
       main = titles[i])
  plot(land,add=TRUE); box()
  plot(filter(clus_ras_vect,cluster==i)$geom, col=col, border=border,
       main = titles[i], add=TRUE)
}

clus_ras_vect$regime = as.factor(clus_ras_vect$cluster)
levels(clus_ras_vect$regime) = titles

# Output levels and mask
write_sf(clus_ras_vect, 'output/regimes.shp')
write_sf(land, 'output/land-utm16.shp')
write_sf(land_lines, 'output/land_lines-utm16.shp')

#function to get extent easily
coords_to_extent_utm = function(xlim, ylim, crs=st_crs(32616)) {
  e = tibble(lat = ylim[c(1,2,2,1,1)],
             lon = xlim[c(1,1,2,2,1)]) %>% 
    st_as_sf(coords = c('lon', 'lat'), crs=4326) %>%  st_bbox() %>% st_as_sfc() %>%
    st_transform(32616) %>% st_bbox() %>% st_as_sfc() %>% st_bbox()
  return(e)
}

# Get an elevation map to underlay
#library(elevatr)
# dem = rast(get_elev_raster(st_transform(land_lines, 4326), z=4))
# hs = shade(terrain(dem, v='slope', unit='radians'),terrain(dem, v='aspect', unit='radians'))
# dem = project(dem, 'epsg:32616')
# dem = crop(dem, extend(clus_ras, 500))
# hs = project(hs, 'epsg:32616')
# hs = crop(hs, extend(clus_ras, 500))
#writeRaster(hs, 'input/hillshade_hires.tif', overwrite=TRUE)
#writeRaster(dem, 'input/dem_hires.tif', overwrite=TRUE)
 
dem = rast('input/dem_hires.tif')
hs = rast('input/hillshade_hires.tif')

#e = coords_to_extent_utm(c(-110,-60), c(10,50))
e = st_bbox(clus_ras_vect)
cluster_map_full = ggplot() + 
  geom_spatraster(data=dem, show.legend = FALSE) +
  scale_fill_gradient(low=grey(0.4), high='#FFFFFF') +
  geom_sf(data=world, size=0.1, fill=NA, color='white')+
  new_scale("fill") +
  geom_sf(data=clus_ras_vect, aes(fill=regime), color=NA) +
  scale_fill_brewer(type='qual', palette=2)+
  geom_sf(data=world, fill=NA, size=0.1, color=grey(0.2)) + 
  coord_sf(xlim=e[c(1,3)], ylim = e[c(2,4)]) +
  theme(legend.position = c(0.85,0.6)) +
  theme(legend.title = element_text( size=6), legend.text=element_text(size=6)) +
  labs(fill='Hurricane regime')

tmp = cluster_map_full + theme(axis.text=element_blank(), axis.ticks = element_blank())
ggsave('figs/cluster_map_full.png', tmp, height=6.5, width=10, dpi=600)  
print(tmp)

#MUlti-panel cluster map

#carribean
e = coords_to_extent_utm(c(-86,-61), ylim=c(14,25))
carrib = cluster_map_full +
  coord_sf(xlim=e[c(1,3)], ylim = e[c(2,4)]) + theme(legend.position='none') +
  ggtitle('C. Caribbean')

#atlantic
e = coords_to_extent_utm(c(-82,-73), ylim=c(25,42))
atlantic = cluster_map_full +
  coord_sf(xlim=e[c(1,3)], ylim = e[c(2,4)]) + theme(legend.position=c(0.75,0.25)) +
  ggtitle('B. Atlantic')

#gulf
e = coords_to_extent_utm(c(-80,-98), ylim=c(18,30))
gulf = cluster_map_full +
  coord_sf(xlim=e[c(1,3)], ylim = e[c(2,4)]) +
  theme(legend.position='none') + 
  ggtitle('A. Gulf')

CLUSTER_FIG = ggarrange(ggarrange(gulf, carrib, nrow=2, heights=c(1.3,1)), atlantic, ncol=2, widths=c(1.3,1))  +
  bgcolor('white')
ggsave('figs/cluster_map_panels.png', CLUSTER_FIG, height=6, width=8, dpi=600)
print(CLUSTER_FIG)


# PLOT 30 YEAR RISK MAPS
riskmap_plot_list = list()
e = st_bbox(clus_ras_vect)
#e = coords_to_extent_utm(c(-80,-98), ylim=c(18,30))
for(ss in 0:5) {
  letter = letters[ss+1] %>% toupper
  label  = bquote(.(letter)*'. '* P[30] ~ Category ~ .(ss))
  if(ss == 0) label  = bquote(.(letter)*'. '* P[30] *' TS-force')
  dat_raster = mask(riskmap_30[[ss+1]], clus_ras_vect)
  if(ss>=2) dat_raster[dat_raster<0.0001] = NA
  
  tmp = 
    ggplot() +
    geom_spatraster(data=dem, show.legend = FALSE) +
    scale_fill_gradient(low=grey(0.4), high='#FFFFFF') +
    geom_sf(data=world, size=0.1, fill=NA, color='white')+
    new_scale("fill") +
    geom_spatraster(data = dat_raster) +
    scale_fill_viridis_c(na.value = NA, option='magma') + 
    theme(legend.position = c(0.85,0.6),
          axis.text=element_blank(), axis.ticks = element_blank()) +
    geom_sf(data=world, fill=NA, size=0.1, color=grey(0.2)) + 
    coord_sf(xlim=e[c(1,3)], ylim = e[c(2,4)]) +
    labs(title=label) +
    labs(fill=expression(P[30]))
  
  riskmap_plot_list[[ss+1]] = tmp
}

#riskfig_2panel = ggpubr::ggarrange(plotlist = riskmap_plot_list[1:2], ncol=1) + bgcolor('white')  
#ggsave('figs/30yr-risk-fig_2panel.png', riskfig_2panel, width=4.5, height=7.6, dpi=600)

riskfig_6panel = ggpubr::ggarrange(plotlist = riskmap_plot_list[1:6], ncol=3, nrow=2) + bgcolor('white')  
ggsave('figs/30yr-risk-fig_6panel.png', riskfig_6panel, width=10.6, height=6, dpi=600)

# PLOT CLUSTER CHARACTERISTICS
risk_summary$cluster = as.numeric(as.character(risk_summary$cluster))
x = clus_ras_vect %>% st_drop_geometry() %>% group_by(cluster) %>% summarize(regime = regime[1])
risk_summary = risk_summary %>% left_join(x)
risk_summary$ss = risk_summary$ss %>% as.factor
levels(risk_summary$ss) = c('TS-force', 1:5)

a = risk_summary %>% ggplot(aes(x=regime, group=ss, fill=regime, y=risk)) +
  geom_bar(aes(alpha=ss), stat='identity', position='dodge', color=grey(0.2)) +
  scale_fill_brewer(type='qual', palette=2) + 
  theme_bw() + 
  theme(legend.position='none') + 
  labs(x=element_blank(), y='30-yr occurrence probability')

t_shift <- scales::trans_new("shift",
                             transform = function(x) {x+8},
                             inverse = function(x) {x-8})
b = risk_summary %>% ggplot(aes(x=regime, group=ss, fill=regime, y=log10(risk+1e-100))) + 
  geom_bar(aes(alpha=ss), stat='identity', position='dodge', color=grey(0.2)) +
  scale_fill_brewer(type='qual', palette=2) + 
  scale_y_continuous(trans = t_shift, breaks = seq(-7,0,1), limits=c(-8,0), labels = paste0('1e',seq(-7,0,1))) +
  labs(x=element_blank(), y=element_blank()) +
  theme_bw()

FIG_CLUSTER_CHARS = ggarrange(a,b, labels=c('A', 'B'), widths=c(0.4,0.6))  
print(FIG_CLUSTER_CHARS)
ggsave('figs/cluster-characteristics.png', FIG_CLUSTER_CHARS, width=8, height=3.5, dpi=600)

#probably move all this to another areas.
#OUTPUT TABLE OF cluster summaries

### Make an ace ggridges plot for each area
# Get area of regimes
regime_area = clus_ras_vect %>% group_by(regime) %>%
  summarize(area_km2=sum(area_km2)) %>% st_drop_geometry()

# Calculate ACE based on Bell 2000
# remove those not on a 6hr basis
# remove those <64kt.
#only include since 1950...
trks = read_csv('input/AL-NECP-combined.csv') %>% 
  filter(max_speed >= 34) %>%
  filter(time %% 600 == 0) %>%
  st_as_sf(coords=c('lon', 'lat'), crs=4326) %>%
  st_transform(st_crs(32616)) %>%
  st_intersection(clus_ras_vect) %>%
  st_drop_geometry() %>%
  filter(substr(date,1,4)>=1949  & substr(date,1,4) < 2022)

date_range = c('1949-01-01 GMT', '2021-12-31 GMT') %>%
  as.Date()
date_df = expand.grid(date=seq(date_range[1],date_range[2], by=1) %>%
                   format('%Y%m%d') %>% as.factor(),
                   max_speed = 0,
                   regime=unique(trks$regime)) %>% as_tibble()

trks = trks %>% select(date, max_speed, regime) %>%
  rbind(date_df)

ace = trks %>% 
  mutate(
    ace_contrib = max_speed ^2,
    year = substr(date,1,4) %>% as.numeric() %>% as.factor(),
    month = substr(date,5,6) %>% as.numeric(),
    day = substr(date,7,8) %>% as.numeric(),
    date = as.Date(paste(year,month,day,sep='-')),
    hurr = max_speed>=64,
    hurr3 = max_speed>=96,
    yday = yday(date),
    week = week(date)) %>%
  arrange(date)

ace_summary = ace %>% group_by(regime, year) %>%
  summarize(ace = sum(ace_contrib),
            hurr_days = sum(hurr),
            hurr3_days= sum(hurr3)) %>% complete(year) %>%
  left_join(regime_area) %>%
  mutate(ace_km2 = ace/area_km2)

regime_summary = ace_summary %>% group_by(regime) %>%
  summarize(Area_km2 = area_km2[1],
            ACE_mn = mean(ace)/Area_km2*1000,
            ACE_sd = sd(ace)/Area_km2*1000,
            hurday_mn = mean(hurr_days)/Area_km2*1000,
            hurday_sd = sd(hurr_days)/Area_km2*1000,
            hur3day_mn = mean(hurr3_days)/Area_km2*1000,
            hur3day_sd = sd(hurr3_days)/Area_km2*1000)
knitr::kable(regime_summary, digits = c(0,0,6,6,6,6,6,6))  
write_csv(regime_summary, 'tabs/hurr-regime-summary.csv')

# Summarize p1, p10, and p30 within clusters
out = list()
out = do.call(rbind, lapply(1:nlyr(riskmap_01), function(i) {
  tmp_out = zonal(riskmap_01[[i]], clus_ras, fun =mean)
  tmp_out$ss = i - 1
  colnames(tmp_out) = c('cluster', 'p1', 'ss')
  tmp_out$regime=titles
  tmp_out = as_tibble(tmp_out)
  tmp_out = tmp_out %>% select(regime, cluster, ss, p1)
  return(tmp_out)
}))
risk_summary = out
risk_summary = risk_summary %>% mutate(
  p10 = 1-(1-p1)^10,
  p30 = 1-(1-p1)^30,
  return_period = 1/p1
) %>% arrange(cluster,ss) %>% select(!cluster)

write_csv(risk_summary, 'tabs/hurr-risk-summary.csv')
knitr::kable(risk_summary, digits = c(0,0,6,6,6,0))  


