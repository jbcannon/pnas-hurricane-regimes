setwd('/home/jcannon/hurrecon-rev')
# Some figures
library(hurrecon)
library(tidyverse)
library(terra)
library(ggplot2)
library(sp)
library(ggpubr)
library(ggnewscale)
data('geographic')

# Load hurricane Michael at landfall
michael = hurrecon::load_hurdat_track('input/AL-hurdat2-1851-2021.csv', 'AL142018')
michael = bind_cols(michael, get_Ha_Hv(michael))
mike_samp = michael %>% slice(c(14,15,16,18:22)) %>%
  mutate(cov = c('W', 'W', 'W', 'L',rep('L',4))) 

landfall = michael[15,]
profiles = get_wind_profiles(landfall)
Vs = hurrecon(mike_samp, 500, 200)
Vs_full = Vs
Vs[Vs<18] = NA

Vmax = c(profiles$summary$Vm_ms[1])
cols = paste0(c('#000000', palette()[2:4]), '99')

#png('figs/hurrecon-wind-profile-demo.png', width=8, height=4, res=600, units='in')
par(mfrow=1:2, mar=c(4,4,1,1))

x = profiles$summary
profile_fit = expand.grid(dir = c('ne','nw','se','sw'),
            R = c(1:300)) %>%
  as_tibble %>% 
  mutate(Rm = x$Rm_km[match(dir, x$dir)], 
         B = x$B[match(dir, x$dir)],
         Vmax = x$Vm_ms[match(dir, x$dir)], 
         Vs = wind_profile_fxn(R, Rm, B) * Vmax)

PROFILE_FIG = profiles$profiles %>% ggplot(aes(y=V * Vmax, x=R)) +
  geom_jitter(aes(color=dir), height=0.00001) +
  theme_classic() +
  geom_line(data = profile_fit, aes(y=Vs, x=R, color=dir), size=1.2, alpha=0.5) +
  labs(x = 'hurricane radius (km)',
       y = bquote(V[s]*', sustained wind speed ('*m~s^{-1}*')')) +
  theme(legend.position = c(0.8,0.8)) +
  scale_color_discrete(name = "direction")
print(PROFILE_FIG)

hurr = raster::raster(Vs) %>% as('SpatialPixelsDataFrame') %>% 
  spTransform('+init=epsg:4326') %>% as.data.frame()
colnames(hurr) <- c("value", "x", "y")

mike_lines = michael %>% group_by(track_name) %>%
  summarize(do_union = FALSE) %>%
  sf::st_cast('LINESTRING')

MICHAEL_VS_FIG = land_lines %>% #filter(name == 'FL') %>%
  ggplot() +
  geom_sf() +
  ggnewscale::new_scale("fill") +
  geom_tile(data = hurr, aes(x=x,y=y,fill=value), height=0.02, width=0.02, alpha=0.8) +
#  geom_sf(data=mike_lines, alpha=0.5, color='white', size=1) +
#  geom_sf(data=mike_lines, alpha=0.5, linetype='dashed') +
  geom_sf_text(data=michael[c(14:16,18:20),], aes(label=time), nudge_x=-1, size=2,color=c(rep('white',4),rep('black',2))) +
  labs(x=element_blank(), y=element_blank(), fill=bquote('wind\nspeed'~(m~s^-1))) +
  theme_classic() +
  theme(legend.position=c(0.25,0.80), legend.background = element_blank()) +
  scale_fill_viridis_c(option='B') +
  theme(legend.key.size=unit(8, 'pt'), legend.title=element_text(size=8)) +
  coord_sf(crs = sf::st_crs(4326),
         xlim=c(-89,-83), ylim=c(26,33)) +
  scale_x_continuous(breaks=c(seq(-88,-84,2)))

plot(MICHAEL_VS_FIG)

HURRECON_DEMO_FIG = ggarrange(PROFILE_FIG, MICHAEL_VS_FIG, labels=c('A', 'B'), widths=c(1.2,1))  +
  bgcolor('white')
print(HURRECON_DEMO_FIG)

ggsave('figs/hurrecon-wind-profile-demo_R2.png',HURRECON_DEMO_FIG,
       width=6, height=3.5, dpi=600)







