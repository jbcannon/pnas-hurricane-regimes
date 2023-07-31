library(tidyverse)
library(sf)
library(terra)
devtools::install_github('jbcannon/hurrecon')
library(hurrecon)

# download hurricane data if needed
hurdat = 'input/AL-NECP-combined.csv'
if(!file.exists(hurdat)) {
  fetch_best_tracks_data(
    path = 'input/AL-hurdat2-1851-2021.csv',
    src='https://www.nhc.noaa.gov/data/hurdat/hurdat2-1851-2021-100522.txt')
  fetch_best_tracks_data(
    path = 'input/NECP-hurdat2-1949-2021.csv',
    src='https://www.nhc.noaa.gov/data/hurdat/hurdat2-nepac-1949-2022-04042023.txt')
  read_csv('input/AL-hurdat2-1851-2021.csv') %>%
    rbind(read_csv('input/NECP-hurdat2-1949-2021.csv')) %>%
    write_csv('input/AL-NECP-combined.csv')
}

# Convert all tracks to shapefile
all_tracks = read_csv(hurdat) %>%
  st_as_sf(coords = c('lon', 'lat'), crs = 4326) %>%
  mutate(ss = cut(max_speed, c(0,34,64,82,96,112,137,500), labels = -1:5),
         basin = substr(track_id,0,2) %>% as.factor())

# Get list of all tracks that reached hurricane status
hurr_list = all_tracks %>% 
  st_drop_geometry %>%
  group_by(track_id) %>%
  summarise(hurr = any(status=='HU')) %>%
  filter(hurr == TRUE) %>%
  pull(track_id)

all_tracks %>% filter(basin == 'EP') %>% pull(date) %>% range()
all_tracks %>% filter(basin == 'AL') %>% pull(date) %>% range()

# Make LINESTRING of all hurricanes
track_lines = all_tracks %>%
  filter(track_id %in% hurr_list) %>%
  group_by(track_id) %>%
  summarize(n = length(status),
            maxSS = max(as.numeric(ss)-1),
            max_speed = max(max_speed),
            do_union=FALSE) %>%
  st_cast('LINESTRING') %>%
  filter(n >1)

# Load land layers
roi = read_sf('input/roi.shp')
land = read_sf('input/land.shp')
land_lines = read_sf('input/land_lines.shp')
land_buff_200km = read_sf('input/land_buff200km.shp')

ggplot() + 
  geom_sf(data=land_buff_200km, linetype='dashed', size=0.1, fill = 'black', alpha=0.1) +
  geom_sf(data=land_lines, mapping = aes(fill=n_10m_d), size=0.1) +
  theme_bw() + theme(legend.position='none') +
  geom_sf(data=roi, linetype='dashed', fill=NA) +
  lims(x=st_bbox(roi)[c(1,3)], y = st_bbox(roi)[c(2,4)]) +
  scale_fill_viridis_d(option='D') +
  theme(panel.background = element_rect(fill = 'lightblue')) +
  geom_sf(data=track_lines, size=0.1)
ggsave('figs/roi-tracks.png', width=5, height=4, dpi=600)

# How many TS/HU make landfall each year AMONG THESE TWO BASINS, by category
track_info = 
all_tracks %>% st_drop_geometry %>% 
  mutate(ss = as.numeric(as.character(ss))) %>%
  mutate(year = as.numeric(substr(date, 0,4))) %>%
  filter(year >= 1950) %>%
  group_by(track_id, basin) %>%
  summarize(mxss = max(ss),
            lf = any(record=='L')) %>%
  mutate(year = as.numeric(substr(track_id,5,8)),
        lf = ifelse(is.na(lf),FALSE, TRUE))

x =track_info %>% 
  group_by(lf, basin, year, mxss) %>%
  summarize(total = length(mxss),
            landfall = sum(lf)) %>%
# make sure all years and ss combinations are represented to factor in zeros
  rbind(expand.grid(lf=c(TRUE,FALSE), 
                    basin=as.factor(c('AL', 'EP')),
                    year=1950:2020, 
                    mxss=0:5,
                    total=0,
                    landfall=0)) %>% 
  group_by(lf, basin, year, mxss) %>%
  summarize(total = sum(total),
            landfall= sum(landfall)) %>%
  group_by(lf, basin, mxss) %>%
  summarize(annual = mean(total)) %>%
  filter(mxss >= 0) %>%
  mutate(lf = as.character(lf))

lf.lbs = c(ALL = 'All TS', `TRUE` = 'Landfall')
data_summary = x %>% group_by(basin, mxss) %>%
  summarize(annual = sum(annual)) %>%
  mutate(lf = 'ALL') %>%
  select(lf, basin:annual) %>%
  bind_rows(filter(x, lf=='TRUE'))

data_summary %>% ggplot(aes(y=annual, x=mxss, fill=basin)) +
  geom_col(position='dodge', color=grey(0.3)) +
    facet_wrap(~lf,
               labeller = labeller(lf = lf.lbs),
               scales='free_y') +
  labs(y='Mean annual tropical cyclones',
       x = 'Intensity (Saffir-Simpson Scale)',
       title='Tropical cyclone frequency',
       subtitle = 'of the Atlantic (AL) and Eastern Pacific (EP) basins (1950-2021)') +
  scale_x_continuous(breaks = c(0:5)) +
  theme_classic() + 
  theme(axis.text.x = element_text(size=12),
        legend.position=c(0.9,0.7)) +
  scale_fill_brewer(type='qual', palette = 3)

ggsave('figs/basin-activity.png', width=6, height=3, dpi=600)

# How many hurricanes by SS category before and after 2004?
all_tracks %>% st_drop_geometry %>%
  filter(track_id %in% hurr_list) %>%
  group_by(track_id) %>%
  summarize(maxSS = max(ss%>%as.character%>%as.numeric),
         year = str_sub(date,1,4)[1],
         period = ifelse(year<2004, 'pre-2004', 'post-2004')[1]) %>%
  count(maxSS,period) %>% pivot_wider(names_from=period, values_from = n)
  
# add Lat/Lon to tracks for modeling wind speed
all_tracks[, c('lon', 'lat')] = all_tracks %>% st_transform(4326) %>% st_coordinates()

# Build a model for each predicting vr at each wind speed

# Helper function to calc wind speed and deal with -999s and 0s
mean_ws = function(x) { # If all -999, then return NA, else, count -999 as 0.
  x = as.vector(x)
  if(all(x<=0)) return(NA)
  x[x<=0] = 0
  return(mean(x))
  }

# Get average wind speed where valid
for(ws in c(34,50,64)) {
  x = all_tracks %>% st_drop_geometry() %>%
    dplyr::select(contains(paste0(ws,'kt'))) %>%
    apply(1, mean_ws)
  var = paste('mean', ws, sep='')
  all_tracks[var] = x
}

all_tracks$year = substr(all_tracks$date,1,4) %>% as.numeric()


#How many observations were used in the shape modeling?
tmp_df = all_tracks %>% filter(track_id %in% hurr_list,
                      mean34 > 0, max_speed >0)
nrow(tmp_df) # n observations?
tmp_df$track_id %>% unique() %>% length()  # n tracks

# how many historical tracks had pressure observatiosn?
all_tracks %>% st_drop_geometry() %>%
  mutate(press_obs = min_press > 0,
         recent = year >= 2004) %>%
  group_by(recent) %>%
  summarize(press= sum(press_obs), n = length(year), p = press/n) 
   
## BUILD MODELS OF WIND SPEED (34, 50, 64) using valid data
## Using speed, pressure, latitude and all 2- and 3-way interactions
mods = list(mean34 = 
              list(p = lm(mean34 ~ (max_speed) * min_press * lat, all_tracks %>%
                            filter(track_id %in% hurr_list,
                                   mean34 > 0,
                                   max_speed > 0,
                                   min_press > 0)),
                   np = lm(mean34 ~ (max_speed) * lat, all_tracks %>%
                             filter(track_id %in% hurr_list,
                                    mean34 > 0,
                                    max_speed > 0))),
            mean50 = 
              list(p = lm(mean50 ~ (max_speed) * min_press * lat, all_tracks %>%
                            filter(track_id %in% hurr_list,
                                   mean50 > 0,
                                   max_speed > 0,
                                   min_press > 0)),
                   np = lm(mean50 ~ (max_speed) * lat, all_tracks %>%
                             filter(track_id %in% hurr_list,
                                    mean50 > 0,
                                    max_speed > 0))),
            mean64 = 
              list(p = lm(mean64 ~ (max_speed) * min_press * lat, all_tracks %>%
                            filter(track_id %in% hurr_list,
                                   mean64 > 0,
                                   max_speed > 0,
                                   min_press > 0)),
                   np = lm(mean64 ~ (max_speed) * lat, all_tracks %>%
                             filter(track_id %in% hurr_list,
                                    mean64 > 0,
                                    max_speed > 0)))
)

# Apply all models
all_tracks$pred34a =ifelse(all_tracks$mean34 > 0 & all_tracks$min_press > 0, predict(mods[[1]][[1]], all_tracks), NA)
all_tracks$pred34b =ifelse(all_tracks$mean34 > 0, predict(mods[[1]][[2]], all_tracks), NA)
all_tracks$pred50a =ifelse(all_tracks$mean50 > 0 & all_tracks$min_press > 0, predict(mods[[2]][[1]], all_tracks), NA)
all_tracks$pred50b =ifelse(all_tracks$mean50 > 0, predict(mods[[2]][[2]], all_tracks), NA)
all_tracks$pred64a =ifelse(all_tracks$mean64 > 0 & all_tracks$min_press > 0, predict(mods[[3]][[1]], all_tracks), NA)
all_tracks$pred64b =ifelse(all_tracks$mean64 > 0, predict(mods[[3]][[2]], all_tracks), NA)

reorg = all_tracks %>% st_drop_geometry() %>%
  filter(year >= 2004) %>%
  filter(max_speed >= 34) %>% 
  dplyr::select(contains(c('pred', 'min_press', 'mean', 'max_speed', 'lat'))) %>%
  mutate(p = ifelse(min_press==-999,'np', 'p'))
df = list()
for(i in 1:nrow(reorg)) {
  row = reorg[i,]
  tmp = tibble(speed = c(34,50,64,34,50,64),
         press = c(rep('p',3), rep('np', 3)),
         predicted = row[,1:6] %>% as.numeric,
         measured = rep(row[,8:10] %>% as.numeric,2),
         max_speed = row$max_speed)
  tmp = tmp %>% filter(!is.na(predicted) & !is.na(measured))
  if(i %% 100==0) cat(i, 'of', nrow(reorg), '\n')
  df[[length(df)+1]] = tmp
}; df = do.call(rbind, df)
df$press = factor(df$press, levels=c('np', 'p'), labels = c('no pressure', 'pressure'))
df$speed = factor(df$speed, levels=c('34','50','64'), labels=paste0('v = ',c('34','50','64')))

df2 = df
data_gap_fill_fig = df2 %>% ggplot(aes(x=predicted, y=measured)) +
  geom_point(alpha=0.1, size=0.4) +
  geom_abline(slope=1, intercept=0, size=0.5, color=grey(0.5), linetype='dashed') +
  facet_grid(rows = vars(press), cols=vars(speed)) +
  theme_classic() +
  theme(panel.border=element_rect(fill=NA)) +
  geom_smooth(method='lm', se=FALSE) +
  labs(x = 'Predicted radius (NM)', y = 'Measured radius (NM)') 
print(data_gap_fill_fig)

mod = mods[[1]]$np
get_mod_info = function(mod){
  tmp = coef(mod)
  tmp = as_tibble(t(tmp))
  names(tmp)[1]='intercept'
  tmp$r2 = summary(mod)$r.squared
  # Get RMSE and n
  rss = c(crossprod(mod$residuals))
  mse = rss/ length(mod$residuals)
  rmse = sqrt(mse)
  tmp$rmse=rmse
  tmp$n = length(mod$residuals)
  # Get P value
  mod_summary = summary(mod)
  tmp$`p-val` = pf(mod_summary$fstatistic[1],              # Applying pf() function
     mod_summary$fstatistic[2],
     mod_summary$fstatistic[3],
     lower.tail = FALSE)
  return(tmp)
}

# Function to make do.call work with pipe
compile = function(list, fn) {return(do.call(fn ,list))}

# Function to match column names in a list of tables
match_names = function(mods) {
  vars = lapply(mods, function(x) colnames(x)) %>% unlist() %>% unique()
  mods = lapply(mods, function(x) {
    for(i in vars) if(!i %in% colnames(x)) x[,i] = NA
    return(x)
  })
  vars = colnames(mods[[1]])
  mods = lapply(mods, function(x) {x = x[,vars]})
  return(mods)
}

radius_models = mods
#output radius model info into a table
TABLE_RADIUS_MODS_PARS = lapply(names(radius_models), function (i) {
  mods = radius_models[[i]]
  mods = lapply(c('p', 'np'), function(x) {
    get_mod_info(mods[[x]]) %>% add_column(var = i, mod=x, .before='intercept')
  })
  mods = match_names(mods)
}) %>% flatten() %>% compile(rbind)
knitr::kable(TABLE_RADIUS_MODS_PARS, digits = c(NA, NA, 1,2,2,2,4,3,3,3,3,1,0,5))


# Estimate hurricane assymetry
asymm = all_tracks %>% filter(status=='HU', year >= 2004) %>%
  select(contains('kt')| contains('mean')) %>% st_drop_geometry()
asymm = lapply(c(34,50,64), function(v) {
  x = asymm %>% dplyr::select(contains(as.character(v))) %>%
    mutate(speed = paste0('v = ', v))
  colnames(x) = c('ne', 'se', 'sw', 'nw', 'mean', 'speed')
  return(x)}) %>% compile(rbind)
x = lapply(1:nrow(asymm), function(x) as.numeric(asymm[x,1:4]) / as.numeric(asymm[x,5])) %>%
  compile(rbind)
colnames(x) = c('ne', 'se', 'sw', 'nw')
asymm = cbind(asymm[, 'speed'], x) %>%
  pivot_longer(2:5, names_to = 'direction', values_to='ratio') %>% 
  filter(!is.na(ratio))
rm(x)

x = asymm
x$speed = 'combined'
asymm = rbind(asymm, x)
fig_asymmetry = asymm %>% ggplot(aes(x=direction, y=ratio, fill=direction)) +
  stat_summary(fun = 'mean', geom='col') +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = .1) +
  facet_wrap(vars(speed), nrow=1) + #scale_fill_brewer(palette='Paired') +
  theme_classic() + 
  theme(legend.position = 'none') +
  geom_hline(yintercept=1, alpha=0.5, linetype='longdash') +
  labs(x='Direction', y='Relative distance') +
  theme(panel.spacing = unit(30, 'pt'),
        panel.background = element_rect(color='black', fill=NA))
plot(fig_asymmetry)

# how many observations
asymm %>% 
  filter(direction=='ne') %>%
  filter(!is.na(ratio)) %>%
  group_by(speed) %>%
  summarize(n = length(ratio))

tab_asymmetry = asymm %>% filter(speed=='combined') %>%
  group_by(direction) %>%
  summarise(ratio_mn = mean(ratio, na.rm=TRUE),
            ratio_sd = sd(ratio, na.rm=TRUE), 

# ---> Save model objects for inclusion in HURRECON model
save(radius_models, tab_asymmetry, file='input/radius_models.RData')

# ---> Output all figures
ggsave('figs/data-gap-fill.png', data_gap_fill_fig, dpi=600,  width=6,height=3.5)

ggsave('figs/radius_assymetry.png', fig_asymmetry, dpi=600,  width=6,height=3.5)


write_csv(tab_asymmetry, 'tabs/assymetry.csv')
write_csv(TABLE_RADIUS_MODS_PARS, 'tabs/radius_mods_pars.csv')

## --->> Additional summary information for writeup
sink('tabs/other_results.txt')
cat('##---> Number of tracks examined that were HURRICANES\n\n')

cat('## Eastern Pacific 1950-2021\n')
grep('EP', hurr_list) %>% length()
cat('## Northern Atlantic 1851-2021\n')
grep('AL', hurr_list) %>% length()

cat('\n\n##---> Number of tracks pre/post 2004\n')
tibble(`yr_gte_2004` = (hurr_list %>% str_sub(-4L, -1L) %>%
                          as.numeric() >= 2004)) %>%
  summary %>% knitr::kable()

cat('\n\n##---> Number of observations with pressure vs no pressure data\n')
all_tracks %>% st_drop_geometry() %>%
  select(basin, year, min_press) %>%
  mutate(recent = year >= 2004) %>%
  mutate(pressure = ifelse(min_press>0,'pressure', 'no pressure')) %>%
  select(pressure, recent, basin) %>% table %>% knitr::kable()
sink()

##FINISHED exploratory analysis... 
# Re run hurrecon in parallel with new extent 


### At the end how many tracks were run...
x  = list.files('data/', pattern='_Vs.tif')
df = tibble(fn=x) %>% mutate(basin = substr(x,1,2))
table(df$basin)
length(x)
