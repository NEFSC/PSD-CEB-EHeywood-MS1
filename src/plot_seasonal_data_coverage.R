library(ggplot2)
library(dplyr)

setwd('C:/Users/Eleanor.heywood/Documents/seal_telemetry/data/')

# Combined deployment locs is the full location dataset filtered by tag dur cutoffdates
dat = read.csv("./L1/locs/Hg_2019-2023_CombinedDeploymentLocs.csv", stringsAsFactors = F)
dat = dat[,c('id', 'datetime', 'sex', 'tag.model.new', 'cutoffstart', 'cutoffend', 'totdatadays')]
dat$datetime = as.POSIXct(strptime(dat$datetime, format = '%Y-%m-%dT%H:%M:%SZ', tz = 'UTC'))
dat$cutoffstart = as.POSIXct(strptime(dat$cutoffstart, format = '%Y-%m-%dT%H:%M:%SZ', tz = 'UTC'))
dat$cutoffend = as.POSIXct(strptime(dat$cutoffend, format = '%Y-%m-%dT%H:%M:%SZ', tz = 'UTC'))

#View(dat[which(is.na(dat$datetime2)),])

dat$monthday = factor(format(dat$datetime, '%m-%d'))
dat$xlabs = format(dat$datetime, '%b-%d')
dat$year = format(dat$datetime, '%Y')
dat$doy = format(dat$datetime, '%j')

dat$ptt = factor(dat$id)


tab <- data.frame(with(data = dat, expr = table(ptt, monthday)))

tab$HasData = ifelse(tab$Freq >= 1, TRUE, FALSE)

meta = dat[, c('ptt', 'sex', 'tag.model.new', 'cutoffstart', 'cutoffend', 'totdatadays', 'year')]
meta = dplyr::distinct(meta)

tab = left_join(tab, meta, by = "ptt")
tab = tab[tab$totdatadays > 10, ]

# Define sex-based outline ranges
datarangesforplot = tab %>% group_by(ptt, cutoffstart, cutoffend) %>% arrange(monthday) %>%
  reframe(#startMD = monthday[which(HasData)[1]-1], 
    startMD2 = monthday[which(HasData)[1]],
    #endMD = monthday[which(HasData)[length(which(HasData))]+1], 
    endMD2 = monthday[which(HasData)[length(which(HasData))]])

tab = left_join(tab, datarangesforplot, 'ptt')

breaks_14d <- levels(tab$monthday)[seq(1, length(levels(tab$monthday)), by = 7)]

bluemale = '#377eb8'
redfemale = '#e41a1c'
# Define a function to get colors based on 'sex'
meta = meta %>% arrange(totdatadays) %>% filter(totdatadays>10)
a = ifelse(meta$sex == 'F', redfemale, bluemale)

p = ggplot() +
  geom_rect(data = NULL, mapping = aes(xmin = stage("01-11"),
                                       xmax = stage("02-29"),
                                       ymin = -Inf, ymax = Inf), fill = 'cornflowerblue', alpha = 0.4) +
  geom_rect(data = NULL, mapping = aes(xmin = stage("02-29"),
                                       xmax = stage("05-31"),
                                       ymin = -Inf, ymax = Inf), fill = 'springgreen', alpha = 0.4) +
  geom_rect(data = NULL, mapping = aes(xmin = stage("05-31"),
                                       xmax = stage("08-31"),
                                       ymin = -Inf, ymax = Inf), fill = 'yellow', alpha = 0.4) +
  geom_raster(data = tab[tab$HasData, ], 
              mapping = aes(x = monthday, 
                            y = reorder(ptt, totdatadays),
                            fill = year), 
              alpha = 0.95) +
  geom_tile(data = tab[tab$HasData, ], 
            mapping = aes(x = monthday, 
                          y = reorder(ptt, totdatadays)),
                          fill = NA, color = 'black') +
  theme_bw() +
  scale_x_discrete(labels = function(x) ifelse(x %in% breaks_14d, x, "")) +
  scale_fill_manual(values = c('gray15', 'gray35', 'gray55', 'gray75', 'gray95')) +
  scale_color_manual(values = c(redfemale, bluemale)) +
  labs(x = "month-day", y = "PTT") +
  theme(axis.text.x = element_text(angle = 75, hjust = 1, color = 'black'), 
        panel.grid.major = element_blank(), 
        axis.text.y = element_text(face = 'bold', colour = a), legend.position = 'inside', 
        legend.position.inside = c(0.9, 0.3)) 
  

p = p + theme(text = element_text(size = 16))

ggsave(filename = '../manuscripts/heywood_etal_2024/AnimalBiotelemetrySubmission/figures/AdditionalFile1_tagdurplot.png',
       plot = p, device = 'png', dpi = 320, 
       width = 170, height = 170, units = 'mm', scale = 1.5)
