# Plot Supplemental Figs. 4 5 and 6
# Plot main manuscript Fig. 4
# EIH
# 2024-09-20

library(readr)
library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)
library(patchwork)
setwd('C:/Users/Eleanor.heywood/Documents/seal_telemetry/')

###############################################################################################################
##################################### SETUP & DATA LOAD #######################################################
###############################################################################################################
# PLOT DIVE STATS INDEPENDENT OF WHETHER OR NOT WE COULD MODEL THEM
###############################################################################################################
##################################### SUPPLEMENTAL FIGURE 5 ###################################################
###############################################################################################################
# Calculate the mean max depth, duration, IDI daily for each PTT
alldives = read_csv('./data/L1/dive/Hg_2019-2023_BEHDiveRecords_QAQC.csv')

# Subset to pre-construction
d = alldives[which(alldives$start_Dive < as.POSIXct('2023-06-01', tz = 'UTC')),]
# Get negative depths
d$Depth = -d$Depth

hist(d$Depth)

# The number of ptts per sex for dive data
tmp = d[,c('ptt', 'sex')]
tmp = distinct(tmp)
table(tmp$sex)


# define color scheme
bluemale = '#377eb8'
redfemale = '#e41a1c'


# Get some time vars
d$doy = as.numeric(format(d$start_Dive, '%j'))
d = d[which(d$doy <= 243), ]
d$md = factor(format(d$start_Dive, '%m-%d'))

labels = distinct(arrange(data.frame(doy = d$doy, md = d$md), doy))
# get the duplicated due to the leap year (2020)
dupes = which(duplicated(labels$doy))-1
labels = labels[-dupes,]

biweeklylabs = labels[seq(1,233, 14), ]
names(biweeklylabs)[2] = 'xlab'

labels = left_join(labels, biweeklylabs, 'doy')
labels = fill(labels, xlab,.direction = 'down')

d = left_join(d, labels, 'doy')

sampsizes = d %>% group_by(sex, xlab) %>% summarise(nptt = n_distinct(ptt), ndives = n_distinct(SDPairID))
sampsizes$pttlab = paste0(sampsizes$sex, '=', sampsizes$nptt)
sampsizes$divelab = paste0(sampsizes$sex, '=', sampsizes$ndives)
biweeklylabs = left_join(biweeklylabs, sampsizes, 'xlab')

gg1 <- ggplot() +
  geom_point(data = d, mapping = aes(x=doy, y = Depth), 
             alpha = 0.25, color = 'gray85',
             position = position_dodge(width = 1)) +
  
  geom_smooth(data = d, mapping = aes(x=doy, y = Depth, color = sex, fill = sex),
              alpha = 0.5, linewidth = 1.5) +
  scale_color_manual(values = c(redfemale, bluemale), labels = c('Female', 'Male')) +
  scale_fill_manual(values = c(redfemale, bluemale), labels = c('Female', 'Male')) +
  
  ylab('max dive depth (m)') +
  theme_bw() +
  scale_x_continuous(name = 'month-day', 
                     breaks = biweeklylabs$doy,
                     labels = biweeklylabs$xlab, limits = c(10,246)) +
  scale_y_continuous(breaks = c(seq(-125, -25, 25), -5)) +
  theme(legend.position = 'inside', 
        legend.position.inside = c(0.065, 0.135),
        legend.title = element_blank(),
        legend.background = element_rect(fill = NA, colour = 'gray20')) +
  coord_cartesian(expand = F, ylim = c(-125, 0))


gg1

gg2 <- ggplot() +
  
  geom_point(data = d, mapping = aes(x=doy, y = DiveDur), 
             alpha = 0.25, color = 'gray85',
             position = position_dodge(width = 1)) +
  
  geom_smooth(data = d, mapping = aes(x=doy, y = DiveDur, color = sex, fill = sex),
              alpha = 0.5, linewidth = 1.5) +
  scale_color_manual(values = c(redfemale, bluemale), labels = c('Female', 'Male')) +
  scale_fill_manual(values = c(redfemale, bluemale), labels = c('Female', 'Male')) +
  
  ylab('dive duration (s)') +
  
  geom_text(data = biweeklylabs, mapping = aes(x = doy, y = rep(c(38,50),17)),
            label = biweeklylabs$pttlab,
            color = rep('black', 34),
            group = 'B',
            size = 2.5,
            hjust = 0,
            nudge_x = 2) +
  geom_text(data = biweeklylabs, mapping = aes(x = doy, y = rep(c(8,20),17)),
            label = biweeklylabs$divelab,
            color = rep('gray50', 34),
            fontface = 'italic',
            size = 2.5,
            group = 'B',
            hjust = 0,
            nudge_x = 1) +
  
  
  theme_bw() +
  theme(legend.position = 'none') +
  scale_x_continuous(name = 'month-day', 
                     breaks = biweeklylabs$doy,
                     labels = biweeklylabs$xlab, limits = c(10,246)) +
  scale_y_continuous(breaks = seq(60, 400, 30)) +

  coord_cartesian(expand = F, ylim = c(0, 400))

gg2
gg3 <- ggplot() +
  geom_point(data = d[d$IDI<=900,], mapping = aes(x=doy, y = IDI), 
             alpha = 0.25, color = 'gray85',
             position = position_dodge(width = 1)) +
  
  geom_smooth(data = d[d$IDI<=900,], mapping = aes(x=doy, y = IDI, color = sex, fill = sex),
              alpha = 0.5, linewidth = 1.5) +
  scale_color_manual(values = c(redfemale, bluemale), labels = c('Female', 'Male')) +
  scale_fill_manual(values = c(redfemale, bluemale), labels = c('Female', 'Male')) +
  
  ylab('inter-dive interval (s)') +
  theme_bw() +
  scale_x_continuous(name = 'month-day', 
                     breaks = biweeklylabs$doy,
                     labels = biweeklylabs$xlab, limits = c(10,246)) +
  scale_y_continuous(breaks = seq(0, 400, 30)) +
  theme(legend.position = 'none') +
  
  coord_cartesian(ylim = c(0,400), expand = F)
gg3

# Combine using patchwork
comb = (gg1 + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), text = element_text(size = 12))) /
  (gg2 + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), text = element_text(size = 12))) /
  (gg3 + theme(axis.title.x = element_blank(), text = element_text(size = 12))) + 
  plot_annotation(tag_levels = 'a')

comb 
ggsave(filename = "./plots/manuscript/FigureS5_divestats_gams.png", plot = comb, 
       device = 'png', dpi = 300, width = 170, height = 170, units = 'mm',scale = 1.5)

###############################################################################################################
###############################################################################################################
###############################################################################################################
