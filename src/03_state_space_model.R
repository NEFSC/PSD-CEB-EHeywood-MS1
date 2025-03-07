# state-space model all valid trips
# EIH
# updated: 2024-09-27-2024
# inputs: HG_2019-2023_HAULOUT_TRIPID_ASSIGNMENT.csv
## libraries
library(readr)
library(magrittr)
library(dplyr)
library(ggplot2)

library(aniMotum)
library(evaluate)
library(parallel)

# use patchwork model fit plot diagnostics
library(patchwork)

rm(list = ls())

################################################## load data ######################################################
setwd('~/seal_telemetry/')    
# 1. load output of 01_preprocess_haulouts_trips.R, sort by id and datetime
# 2. get qualifying trips (bookended by a confirmed haul out bout)
# 3. calculate segment length and subset segments with >= 50 positions for model convergence
# 4. fit state-space model with a continuous-time correlated random walk movement process in aniMotum

# load data
df = read_csv("./data/L1/locs/HG_2019-2023_HAULOUT_TRIPID_ASSIGNMENT.csv") %>%
  arrange(id, datetime)

# look at the location ping rate 
pingrate = df %>% filter(!is.na(SegID)) %>% group_by(SegID) %>% arrange(datetime) %>%
  mutate(tdiff = difftime(time1 = datetime, time2 = lag(datetime, n = 1), tz = "UTC", units = 'sec')) %>% 
  dplyr::select(SegID, datetime, tdiff) %>% ungroup() %>%
  group_by(SegID) %>%
  summarise(meantdiff = mean(tdiff, na.rm=T))

# get mean and sd ping rate in hours 
print(round(as.numeric(mean(pingrate$meantdiff, na.rm=T))/60/60, digits = 2))
print(round(as.numeric(sd(pingrate$meantdiff, na.rm=T))/60/60, digits = 2))

## get qualifying trips (bookended by confirmed haulout periods)
df$eventID = NA
df$eventID[is.na(df$HauloutID)] = paste("T-", df$TripID[is.na(df$HauloutID)], sep="")
df$eventID[is.na(df$TripID)] = paste("H-", df$HauloutID[is.na(df$TripID)], sep = "")

hoL = split(df, df$id)
qual_trip_ids = c()
for (i in 1:length(hoL)){
  tmp = hoL[[i]]
  tmp = arrange(tmp, datetime)
  events = unique(tmp$eventID)
  ht = sapply(strsplit(events, "-"), "[", 1)
  # get indices of H and T
  ts = which(ht == "T")
  hs = which(ht=="H")
  print(all(diff(ts) == 2)) # ok good, these are all true
  print(all(diff(hs) == 2))
  # if the last trip index is greater than the last haulout index... then it does not qualify
  if (ts[length(ts)] > hs[length(hs)]){
    quals = ts[1:length(ts)-1]
  }else{quals = ts}
  
  qual_ids = events[quals]
  qual_ids = gsub("T-", "", qual_ids)
  qual_trip_ids = c(qual_trip_ids, qual_ids)
  rm(qual_ids, quals, hs, ts, ht, events, tmp)
}

# subset for at-sea locations
df_ = df[!df$final_haulout, ]

# get number of segments per tripqualifying trip 
seg_sum = df_ %>% group_by(TripID) %>% summarise(NSegTotal = length(unique(SegID)))

# get seg length and subset to segments >= 50
df_ = df_ %>% group_by(SegID) %>% 
  mutate(SegLength = n())%>%
  filter(SegLength >= 50)

# get number of segments per trip after eliminating segments < 50 obs
seg_sum2 = df_ %>% group_by(TripID) %>% summarise(NSegQual = length(unique(SegID)))

segsum = left_join(seg_sum2, seg_sum, 'TripID')

# Get only completely modelable trips
completetripids = segsum$TripID[which(segsum$NSegQual == segsum$NSegTotal)]

# downstream: will use qual_trip_ids + completetripids to subset the modeled trips for any TRIP Level analysis that requires complete qualifying trips... 
# for dive analyses, we still want modeled partial trips.

# get number of modeled positions for complete and qualifying trips
completequaltrips = intersect(completetripids, qual_trip_ids)

# Compare the number of positions modeled for completed, qualifying trips in the pre-construction phase analysis
analysespos = df_ %>% filter(TripID %in% completequaltrips) %>% group_by(TripID) %>% mutate(tripend = max(datetime)) %>%
  filter(tripend < as.POSIXct(strptime(x = '2023-06-01 00:00:00', format = '%Y-%m-%d %H:%M:%S', tz = 'UTC')))
n_distinct(analysespos$id)
n_distinct(analysespos$TripID)

# to the total possible positions (those positions from all trips ending prior to the construction date cutoff)
totalpos = df_ %>% group_by(TripID) %>% mutate(tripend = max(datetime)) %>%
  filter(tripend < as.POSIXct(strptime(x = '2023-06-01 00:00:00', format = '%Y-%m-%d %H:%M:%S', tz = 'UTC')))

round(nrow(analysespos) / nrow(totalpos)*100) # 90%
########################################################### fit ssm-crw #####################################################################

# model each trip segment
anilocs = df_[,c('SegID', 'datetime', 'lc','lon', "lat", "error.radius", 'smaj', "smin", "eor")]

# segdurs = anilocs %>% group_by(SegID) %>% summarise(segdur = as.numeric(difftime(max(datetime), min(datetime), tz = 'UTC', units = 'days')))
# segdurs$ptt = sapply(strsplit(segdurs$SegID, '-'), '[[', 1)
# meansegdursperptt = segdurs %>% group_by(ptt) %>% summarise(meansegdur = mean(segdur))

colnames(anilocs)[c(1,2)] = c('id', 'date')
aniList = split(anilocs, anilocs$id)

# clean workspace
rm(seg_sum, seg_sum2, segsum, hoL)

# fit ssm-crw 
set.seed(59)
# set time step for regularization
ts = 2
vmax = 3.5 # m/s Gallon et al. 2007
segfits = fit_ssm(x = anilocs, 
                  vmax = vmax, 
                  model = 'crw', 
                  time.step = ts,
                  control = ssm_control(optim = 'optim', method = 'L-BFGS-B', verbose = 0))

s = as.data.frame(summary(segfits)$Stattab)
any(s$converged == FALSE) # FALSE - everything converged

# look at standard error estimations that were not calculated
idsforrefit = sapply(X = 1:nrow(segfits), FUN = function(x){
  tmppar = as.data.frame(segfits$ssm[[x]]$par)
  
  if (any(is.na(tmppar$`Std. Error`))){
    idforrefit = segfits$id[x]
  }else{idforrefit = NA}
  
  return(idforrefit)
})

idsforrefit = idsforrefit[!is.na(idsforrefit)]

# refit any that failed to converge or where standard error estimates failed
# fit without estimating psi
if (length(idsforrefit) > 0){
  refitanilocs = anilocs[anilocs$id %in% idsforrefit, ]
  segrefits = fit_ssm(x = refitanilocs, 
                      vmax = vmax, 
                      model = 'crw', 
                      time.step = ts, 
                      control = ssm_control(optim = 'optim', method = 'L-BFGS-B', verbose = 0),
                      map = list(psi = factor(NA)))
  
  srefit = summary(segrefits)
  # replace with the refitted now that the standard error issues are resolved
  segfits[which(segfits$id %in% idsforrefit), ] = segrefits
}

# save for later validation & plotting
save(segfits, file = './data/L2/SSM/Hg-2019-2023-SSM-ModelObjects.RData')

########################################## reroute around land ################################
fitrr = route_path(segfits, what = 'predicted', map_scale = 10, buffer = 10000)

# get predictions and ptt id #s from the segID string
preds_reroute = grab(x = fitrr, what = 'rerouted')
preds_reroute$ptt = sapply(strsplit(preds_reroute$id, '-'), '[', 1)

preds = grab(x = segfits, what = 'predicted')
preds$ptt = sapply(strsplit(preds$id, '-'), '[', 1)


### create leaflet dynamic map to check the fits and make sure they are sensical 

library(leaflet)
library(htmlwidgets)
m = leaflet(preds_reroute) %>%
  setView(lng = mean(preds_reroute$lon), lat = mean(preds_reroute$lat), zoom = 5) %>%
  addProviderTiles(providers$Esri.WorldImagery)

preds_reroute$datenum = as.numeric(preds_reroute$date)

for (i in 1:length(unique(preds_reroute$ptt))){
  id = unique(preds_reroute$ptt)[i]
  tmp = preds_reroute[preds_reroute$ptt == id,]
  pal <- colorNumeric(palette = 'viridis', domain = c(min(tmp$datenum), max(tmp$datenum)))
  
  m = m %>% 
    addCircleMarkers(data = tmp,
                     lng = ~lon, 
                     lat = ~lat, 
                     color = ~pal(datenum),
                     label = tmp$date, 
                     stroke = TRUE,
                     fillOpacity = 0.5, 
                     radius = 5,
                     group = as.character(unique(tmp$ptt)))
}

# Add layers control
m = m %>% # Layers control
  
  addLayersControl(overlayGroups = c(unique(preds_reroute$ptt)),
                   position = c("bottomright"),
                   options = layersControlOptions(collapsed = TRUE)
                   )
m

saveWidget(widget = m, file = "./plots/ssm_validation/ssmcrw_reroute_validationmap.html")

# write out prediction dataframes
write_csv(x = preds_reroute, file = "./data/L2/SSM/Hg-2019-2023-SMM-Tracks_reroute.csv")
write_csv(x = preds, file = "./data/L2/SSM/Hg-2019-2023-SMM-Tracks.csv")

############### subset qualifying and complete trips ##########################. 
# qual_trip_ids - are all trips that are bookended by a haulout
# completetripids - are all trips that have complete modelable segments (all segments of trip have >= 50 locations)

# for regular preds
preds$TripID = paste(preds$ptt, sapply(strsplit(preds$id, '-'), '[', 2), sep = "-")
predssub = preds[which(preds$TripID %in% qual_trip_ids), ]
predscomplete = predssub[which(predssub$TripID %in% completetripids), ]

# for rerouted preds
preds_reroute$TripID = paste(preds_reroute$ptt, sapply(strsplit(preds_reroute$id, '-'), '[', 2), sep = "-")
preds_reroutesub = preds_reroute[which(preds_reroute$TripID %in% qual_trip_ids), ]
preds_reroutecomplete = preds_reroutesub[which(preds_reroutesub$TripID %in% completetripids), ]

# write out only completed trips
write_csv(x = predscomplete, file = "./data/L2/SSM/Hg-2019-2023-SMM-Tracks-CompleteTrips.csv")
write_csv(x = preds_reroutecomplete, 
          file = "./data/L2/SSM/Hg-2019-2023-SMM-Tracks-CompleteTrips-Reroute.csv")
# end