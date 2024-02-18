#--- GPS data cleaning prototype script ----------------------------------------#
#--- 2023-01-29 ----------------------------------------------------------------#
#--- Original author: Dani Berger (danielle.berger@usu.edu) --------------------#
#--- Adapted by Kezia Manlove (kezia.manlove@usu.edu) for SCV2 data structure --#

# Script overview: 
# (Each section can be accessed using navigation drop-down located in the slider
# between this pane and the panel below)

# I. Load required libraries and data
# II. Annotate gps data with animal- and capture-specific fields and define dates
# III. Project data to lat-long
# IV. Drop points inside bounding boxes around collar storage locations
# V. Remove fixes before collar deployment or after last good data date
# VI. Reformat the data for AMT
# VII. amt cleaning steps
### A. Remove fixes close to time of capture
### B. Remove fixes with high DOP
### C. Flag duplicates
### D. Flag Fast Steps
### E. Flag Fast Roundtrips
### F. Remove flagged rows
# VIII. Convert the dataset to a df and save as csv
# IX. Build basic maps of seropos's and seroneg's



# I. Load required libraries and data ----
library(devtools)
devtools::install_github("jmsigner/amt") # pull in dev version of amt,
## which has better cleaning functions (esp the fast round trips function)
library(amt)       # dev amt to id fast round trips
library(tidyverse) # tidyverse to script with amt
library(sf)        # sf to 
library(lubridate) # lubridate to manage dates
library(mapview)   # mapview to pull in desired tiles
library(leaflet)   # leaflet to build pretty maps

# Read in datasets
locs_in <- read.csv("ExampleScripts/Data/UDWR_Elk_BookCliffs.csv", header = T)
# checking here to see whether any animal wore two different collars
locs_in$animal_collar <- paste(locs_in$uniqueID, "_", locs_in$collarID, sep = "")
# cutting down to just one row per animal-collar combo
animals_in <- locs_in[!duplicated(locs_in$animal_collar), ]
# there are no extras (nrow(animals_in) is the same if I look for duplicates of
# animal_collar combos or just uniqueIDs)

#locs_in <- read.csv(file = "ExampleScripts/Data/SCV2DeerExample/LocationData/Nebo_thru_feb23.csv", 
#                    header = T)
# animals_in <- read.csv("ExampleScripts/Data/SCV2DeerExample/AnimalData/AnimalData_in.csv", 
#                        header = T)
# deploys_in <- read.csv("ExampleScripts/Data/SCV2DeerExample/CollarDeploymentData/CollarDeploymentData_in.csv", 
#                        header = T)
# caps_in <- read.csv("ExampleScripts/Data/SCV2DeerExample/CaptureData/Utah_2023-01-27 - Utah_2023-01-27.csv", 
#                     header = T)

# II. Annotate gps data with animal- and capture-specific fields and define dates ----
## A. Annotate trajectories with animal-level fields ----
### i. Pull out top row of collar data for each animal
# num_collars <- length(levels(factor(locs_in$collarID)))
# locs_in$collarID_char <- as.character(locs_in$collarID)
# collar_vec <- vector("list", num_collars)
# for(i in 1:num_collars){
#   k <- subset(locs_in, collarID_char == levels(factor(locs_in$collarID))[i])
#   collar_vec[[i]] <- k[1, ]
# }

### ii. Add collar id and corresponding animal data onto deployment data
# KRM specific: drop Deer41 from La Sals
# deploys_in <- subset(deploys_in, animal_id != "UT_LS_SCV2_Deer41")
# # build a deployment ID field containing the animal id and the collar deployment info
# deploys_in$animal_deployment_id <- paste(deploys_in$animal_id, "_", 
#                                          deploys_in$collar_deployment_id, sep = "")
# 
# ### iii. build birth_date field for animal_dat
# ##### I'm assigning all animals born in a year a June 1 birth date
# animals_in$birth_day_char <- paste(floor(animals_in$birth_cohort), 
#                                    "-06-01 12:00:00", sep = "") 
# ### convert that birth day to posixct
# animals_in$birth_day <- lubridate::ymd_hms(animals_in$birth_day_char)
# 
# ### iv. loop over collars, pull animal data on each collar, and append to deployment data
# # build a storage object to write the annotated data into
# deploys_update <- vector("list", num_collars)
# for(i in 1:num_collars){
#   # pick off all deployments for each collar
#   deploys_update[[i]] <- subset(deploys_in, levels(factor(locs_in$collarID))[i] == 
#                                   as.character(collar_id))
#   # only do this part if any deployments are located.
#   if(nrow(deploys_update[[i]]) != 0){ 
#     # pull off the animal ID from the deployment data
#     animal_id_in <- as.character(deploys_update[[i]]$animal_id[1]) 
#     # pull the corresponding animal data based on the animal ID
#     animal_dat <- subset(animals_in, as.character(animal_id) == as.character(animal_id_in))     
#     # append the herd name from animal_dat onto the deployment data
#     deploys_update[[i]]$herd <- rep(animal_dat$herd, nrow(deploys_update[[i]]))
#     # append the sex from animal_dat onto the deployment data
#     deploys_update[[i]]$sex <- rep(animal_dat$sex, nrow(deploys_update[[i]]))    
#     # append the birth cohort from animal_dat onto the deployment data
#     deploys_update[[i]]$birth_day <- rep(animal_dat$birth_day, nrow(deploys_update[[i]]))     
#     # append the area name from animal_dat onto the deployment data
#     deploys_update[[i]]$area <- rep(animal_dat$area, nrow(deploys_update[[i]]))  
#     }
# }
# 
# # bind all submatrices within the deploys_update list together into a dataframe. 
# deploys_update2 <- do.call("rbind", deploys_update)
# 
# # if end_date is blank, fill in with today's date.
# deploys_update2$end_date[which(is.na(deploys_update2$end_date))] <- as.character(print(now()))

# convert locs_in$date to posixct
locs_in$dateYearAndJulian <- lubridate::ymd_hms(as.character(locs_in$dateYearAndJulian))
# convert deployment date and end date to posixct
# deploys_update2$deployment_date <- lubridate::ymd(deploys_update2$deployment_date)
# deploys_update2$end_date <- lubridate::ymd_hms(deploys_update2$end_date)
# deploys_update2$birth_day <- lubridate::ymd_hms(as.character(deploys_update2$birth_day))

# ## B. Annotate trajectories with deployment data ----
# deploys_list <- vector("list", nrow(deploys_update2))
# for(i in 1:nrow(deploys_update2)){
#   # cut down to fixes from the focal collar
#   k_collar <- subset(locs_in, collarID_char == as.character(deploys_update2$collar_id[i]))
#   # cut down to fixes in the focal deployment (within the collar)
#   deploys_list[[i]] <- subset(k_collar, dateYearAndJulian > deploys_update2$deployment_date[i] &
#                                 dateYearAndJulian < deploys_update2$end_date[i])
#   if(nrow(deploys_list[[i]]) >= 1){ # if we have at least one point for this deployment
#     # bind on the following fields with animal-level values replicated across all fixes
#     deploys_list[[i]]$animal_id <- rep(deploys_update2$animal_id[i], nrow(deploys_list[[i]]))
#     deploys_list[[i]]$sex <- rep(deploys_update2$sex[i], nrow(deploys_list[[i]]))
#     deploys_list[[i]]$animal_deployment_id <- rep(deploys_update2$animal_deployment_id[i], 
#                                                   nrow(deploys_list[[i]]))
#     deploys_list[[i]]$age <- deploys_list[[i]]$dateYearAndJulian - deploys_update2$birth_day[i]
#   }
# }
# 
# # bind deploys_list submatrices into a single dataframe
# gps_prepped <- do.call("rbind", deploys_list)
gps_prepped <- locs_in

## C. Annotate trajectories with capture-level data ----
### i. Set up dates in caps_in
### My dates are just date, not time, so I'm setting the time to 12:00:00 for that date
### This allows for comparability across the different POSIXct dates

# paste " 12:00:00" onto the end of the capture date in all cases
# caps_in$Collection.Date <- ifelse(caps_in$Collection.Date == "", NA,
#                                   paste(caps_in$Collection.Date, " 12:00:00", sep = "")) 
# # convert to posixct
# caps_in$Collection.Date <- lubridate::mdy_hms(as.character(caps_in$Collection.Date)) 
# # same thing for birth day
# animals_in$birth_day <- lubridate::ymd_hms(animals_in$birth_day_char) 
# 
# ### ii. Loop over animal IDs in gps_prepped
# animal_ids <- levels(factor(gps_prepped$animal_id))
# animal_gps <- vector("list", length(animal_ids))
# for(i in 1:length(animal_ids)){
#   # pull off all fixes for the focal animal
#   animal_gps[[i]] <- subset(gps_prepped, animal_id == animal_ids[i])
#   # pull off all testing events for the focal animal
#   animal_caps <- subset(caps_in, as.character(Animal.ID) == animal_ids[i])
#   # add empty storage vectors onto the gps dataframe
#   animal_gps[[i]]$serol_cat <- animal_gps[[i]]$serol_quant <- rep(NA, nrow(animal_gps[[i]]))
#   animal_gps[[i]]$days_since_cap <- rep(NA, nrow(animal_gps[[i]]))
#   # loop over capture events to population gps data fields
#   if(nrow(animal_caps) >= 1){ # if we have at least one capture record per animal...
#     for(j in 1:nrow(animal_caps)){ # loop over each capture event within an animal
#       # subset down to all locations on or after the focal capture date
#       cap_day <- which(animal_gps[[i]]$dateYearAndJulian >= animal_caps$Collection.Date[j])
#       # fill in the categorical serological result
#       animal_gps[[i]]$serol_cat[cap_day] <- as.character(animal_caps$Result_sero[j])
#       # fill in the quantitative serological result
#       animal_gps[[i]]$serol_quant[cap_day] <- animal_caps$Pcinh[j]
#       animal_gps[[i]]$days_since_cap[cap_day] <- (animal_gps[[i]]$dateYearAndJulian[cap_day] - 
#                                                     animal_caps$Collection.Date[j])/24
#       # (default unit for time differences is in hours, so "/24" rescales to days)
#     }
#   }
#   print(i)
# }
# 
# gps_prepped2 <- do.call("rbind", animal_gps)
# head(gps_prepped2)
# 
# # rename gps_prepped2 "dat_in"
gps_prepped2 <- gps_prepped
dat_in <- gps_prepped2

# could save dat_in at this point if your run time up to here is slow. 


# III. Project data to lat-long ----
## A. Generic specification of sf object with lat-long ----
# Prototype data are in lat-long, so that's the projection I'm applying here. 
dat_in_sf <- st_as_sf(x = dat_in,
                      coords = c("longitude", "latitude"),
                      crs = 4326)

# 32612 is crs for utms for the La Sals, my example dataset.
# 4326 is crs for lat-long

# Make the point geometry into two separate coordinate columns
dat_in_coords <- as.data.frame(st_coordinates(dat_in_sf))

# could plot these out just to see how they look, but plot construction is
# quite slow due to the number of points. 
# plot(dat_in_coords[, 1] ~ dat_in_coords[, 2])

## B. Transforming data from other crs's ----
## This section is commented out since it does not apply to Utah sites.
## To use, highlight this section and remove one level of comments (use all 
## lines with just one # )
## EXAMPLE: Transform lat/long data to UTM
## crs here is the epsg code for the coordinate reference system your
## data come in (32612 is the crs for UTMs in Utah)
# dat_utmtransform <- st_transform(dat_in_sf, crs = 32612)

## Make the point geometry into two separate coordinate columns
# dat_coords <- as.data.frame(st_coordinates(dat_utmtransform))

## Cbind back to the dat_in_sf
## This step is necessary because the key_field is not actually unique at this
## point. We still have temporal duplicates and will need to join on lat/long as
## well. 
# dat_latlong_all <- cbind(dat_in_sf, dat_coords)

# IV. Drop points inside bounding boxes around collar storage locations ----
## A. set box coordinates around removal zones ----
## This section is commented out in this script because I don't have locations
## to drop for now. Comment out the single-# lines to use. 
## EXAMPLE: 
## Remove points from: 
## 1. Zion research trailer

# ZRT_p1_latlong = st_point(c(37.208576, -112.980949))
# ZRT_p2_latlong = st_point(c(37.208605, -112.980771))
# ZRT_p3_latlong = st_point(c(37.208404, -112.980735))
# ZRT_p4_latlong = st_point(c(37.208377, -112.980901))
# 
## Zion main office (just so we have a second point to remove for demo)
# ZMO_p1_latlong = st_point(c(37.209182, -112.981125))
# ZMO_p2_latlong = st_point(c(37.208882, -112.980558))
# ZMO_p3_latlong = st_point(c(37.20948488235225, -112.98019562476699))
# ZMO_p4_latlong = st_point(c(37.20965551930826, -112.98051549480576))
# 
## B. Create polygons from list of sf points in A. ----
## create polygon around Zion research trailer
# ZRT_poly = st_multipoint(c(ZRT_p1_latlong, ZRT_p2_latlong, ZRT_p3_latlong, ZRT_p4_latlong)) %>%
#   st_cast("POLYGON") %>%
#   st_sfc(crs = 4326) 
# mapview(ZRT_poly)
# 
## create polygon around Zion main office
# ZMO_poly = st_multipoint(c(ZMO_p1_latlong, ZMO_p2_latlong, ZMO_p3_latlong, ZMO_p4_latlong)) %>%
#   st_cast("POLYGON") %>%
#   st_sfc(crs = 4326) 
# mapview(ZMO_poly)
# 
## C. Combine polygons into one sf object ----
## combine polygons
# removal_polys_utms <- c(ZRT_poly, ZMO_poly)
## name polygons
# names(removal_polys_utms) <- c("ZRT_poly", "ZMO_poly")
# 
## D. Project polygons to utms ----
# removal_polys_utms <- st_transform(removal_polys_utms, crs = 32612)
# 
## assess polygons on map
# mapview(removal_polys_utms)
## write to shape file
# st_write(removal_polys_utms, "SampleOutput/Spatial/GPS_data_cleaning/removal_polys_utms.shp", append = FALSE)
# 
## E. Create simple feature data frame of points in bounding box ----
# dat_in_sf <- dat_in %>%
#   st_as_sf(coords = c("EAST","NORTH"), crs = 32612, na.fail = FALSE,
#            remove = FALSE)
# 
## Create sf data frame with points only within bbox
# filtered_removal_polys <- st_filter(dat_in_sf, removal_polys_utms) %>%
#   mutate(within_bbox = "1")  # assign points in bbox a value of 1
# mapview(removal_polys_utms) + filtered_removal_polys
## IF NO POINTS ARE CAUGHT BY THIS FILTER, YOU'LL GET AN ERROR HERE. 
# 
## F. Remove points within bounding box ----
# dat_in_sf <- dat_in_sf %>%
#   # assign values outside of bbox = zero. Points within bbox = 1 from above
#   mutate(out_bbox = ifelse(!X %in% filtered_removal_polys$X, "0", "1")) %>%
#   filter(out_bbox == 0) # remove points within bbox
# 
## Convert back to a regular data frame 
# dat_in <- dat_in_sf %>% 
#   st_drop_geometry() %>% 
#   select(-out_bbox)




# V. Remove fixes before collar deployment or after last good data date ----
## A. Read in deployment data
#deploys_in$dt <- lubridate::ymd(deploys_in$deployment_date)
#deploys_in$deployment_date <- ymd(deploys_in$deployment_date)
#deploys_in$end_date <- ymd(deploys_in$end_date)

## B. Break animals into list and check each deployment in isolation
# dat_list_trunc <- vector("list", length(levels(factor(dat_in$animal_deployment_id))))
# deployments <- levels(factor(dat_in$animal_deployment_id))
# for(i in 1:length(deployments)){
#   deployment_locs <- subset(dat_in, animal_deployment_id == deployments[i])
#   deployment_dat <- subset(deploys_update2, animal_deployment_id == deployments[i])
#   deployment_locs$collar_deployment_date <- rep(deployment_dat$deployment_date, nrow(deployment_locs))
#   deployment_locs$collar_end_date <- rep(deployment_dat$end_date, nrow(deployment_locs))
#   too_early <- which(deployment_locs$dateYearAndJulian <= (deployment_dat$deployment_date + days(2)))
#   too_late <- which(deployment_locs$dateYearAndJulian >= (deployment_dat$end_date - days(1)))
#   removals <- c(too_early, too_late)
#   dat_list_trunc[[i]] <- deployment_locs[-removals, ]
# }
# 
# gps_trimmed <- do.call("rbind", dat_list_trunc)
# # check how many observations are retained
# dim(gps_trimmed)

gps_trimmed <- dat_in


# VI. Reformat the data for AMT ----
# Make GPS data into a track for AMT
data_track <- gps_trimmed %>% 
  make_track(.x = longitude, 
             .y = latitude, 
             .t = dateYearAndJulian, # previously "dt" 
             crs = 4326, # recall, my data are in lat-long
             dop = dop, # my dop field is named "dop"
             # Key_Field = X, 
             # These arguments just append extraneous fields that exist in the 
             # gps_trimmed df onto the appropriate tracks. The names of the 
             # following fields can agree exactly with the corresponding field 
             # name in gps_trimmed.
             animal_id = uniqueID, 
             sex = sex,
             age = currentAge,
             cap_area = realCaptureArea
             # serol_cat = serol_cat,
             # serol_quant = serol_quant,
             # collar_deployment_date = collar_deployment_date, 
             # collar_end_date = collar_end_date
  )

# Group the track by animal_id and arrange chronologically within each individual
data_track <- data_track %>% 
  mutate(animal_id = as.factor(animal_id)) %>%
  group_by(animal_id) %>% 
  arrange(t_, .by_group = TRUE) %>% 
  ungroup()

# amt likes nested data frames, so nest tracks by animal_id
# (the tidy vsn is commented out; the one-liner in 322 appears to get to the same place)
# data_track_nest <- data_track %>% 
#   # These are all the columns that should go in the data field
#   # We will have one row for each Animal_ID and a list of data for each Animal_ID
#   # "nested" in the second column
#   nest(data = c(x_, y_, t_, 
#                 dop, 
#                 animal_id,
#                 collar_deployment_date, 
#                 collar_end_date
#   ))
data_track_nest <- data_track %>% nest(data = -"animal_id")

# check the data for the first individual
data_track_nest$animal_id[[1]]




# VII. amt cleaning steps ----
## A. Remove fixes close to time of capture ----
# Remove all fixes for the first 3 days following capture. 
data_track_nest <- data_track_nest %>%
  mutate(data = map(data, function(x)
    x %>% remove_capture_effect(start = days(3))))

## B. Remove fixes with high DOP ----
# Plot data to see what gets lost by filtering by DOP (to figure out what would be
# a reasonable "high" value)
# Make a sequence of dop values from 1 to max
dop_val <- seq(from = 1, to = max(gps_trimmed$dop, # update dat_in$dop to the name of your DOP field
                                  na.rm = TRUE), by = 1)  

# Make an empty vector the same length as the dop_val sequence
dop_out <- rep(NA, length(dop_val))

# Count the number of fixes with DOP greater than index i
for(i in 1:length(dop_out)){
  dop_out[i] <- sum(dat_in$dop >= i, na.rm = TRUE)
}

# Combine two vectors into a data frame for plotting 
dop_plot <- as.data.frame(cbind(dop_val, dop_out))

# Plot the DOP cutoff vs. the number of rows we would discard if we used that 
# cutoff value. 
ggplot(dop_plot, aes(x = dop_val, y = dop_out)) +
  geom_point() +
  geom_line()

# Here, I use a dop cut-off of 3 and keep all rows that contain NAs. 
data_track_nest <- data_track_nest %>%
  mutate(data = map(data, function(x)
    x %>% filter(dop <= 3 | is.na(dop))))



## C. Flag duplicates ----
# Flag any points within 10 minutes of another fix
# THIS IS A LITTLE SLOWER -- ABOUT 20 SECONDS ON THIS DATASET
data_track_nest1 <- data_track_nest %>%
  mutate(data = map(data, function(x)
    x %>% flag_duplicates(gamma = minutes(10))))

# # R doesn't like that there are NAs in the DOP column
# # Since we filtered out all DOPs over 20, we could fill in NAs with identical 
# # large values (if we are comparing two fixes that are both missing DOP)
# 
# data_track_nest1a <- data_track_nest1 %>% 
#   mutate(data = map(data, function(x)
#     x %>% replace_na(list(dop = 25))))

# check what this flags
# data_track_nest1b <- flag_duplicates(data_track_nest1, gamma = minutes(10),
#                                           time_unit = "mins", DOP = "dop")
# data_track_nest1c <- data_track_nest1a %>%
#   mutate(data = map(data, function(x)
#     x %>% flag_duplicates(gamma = minutes(10))))

# Looking at a few individuals, I'm fairly comfortable saying that DOP wasn't missing
# sporadically, so I'm going to continue to use my work-around by substituting in a
# large dop value when dop is missing. If dop is identical for all the fixes, then
# the function compares the distance between the previous fix and the duplicate 
# fixes and keeps the closest one. 
#p <- data_track_nest1$data[[1]]

# I need to remove the temporal duplicates before proceeding. 
data_track_nest2 <- data_track_nest1 %>%  # update original object before pipe to 
  # data_track_nest1c if you use the commented-out functions above
  mutate(data = map(data, function(x)
    x %>% filter(duplicate_ == FALSE)))



## D. Flag Fast Steps ----
# Distance is assumed to be in kilometers and time must be a lubridate object
delta <- calculate_sdr(speed = 3, # speed given in km/hr
                       time = minutes(60))
# delta <- 2500. These values are totally made up. Basically, I'm saying that a fast-
# moving mule deer would go 3km in an hour. 

# Assuming speed is provided in km/h.
data_track_nest3 <- data_track_nest2 %>%
  mutate(data = map(data, function(x)
    x %>% flag_fast_steps(delta = delta)))





## E. Flag Fast Roundtrips ----
# The recommended default epsilon is 10. 

data_track_nest4 <- data_track_nest3 %>%
  mutate(data = map(data, function(x)
    x %>% flag_roundtrips(delta = delta, 
                          epsilon = 40, 
                          time_unit = "hours")))

# check to see if any points were flagged. 
# data_track_nest5 <- data_track_nest3 %>%
#   mutate(data = map(data, function(x)
#     x %>% filter(fast_roundtrip_ ==TRUE)))
# 
# print(data_track_nest2, n = 405)
# 
# (for me, nothing got flagged here)

## F. Remove rows that have TRUE flags ----
### i. Unnest the track ----
### Remove individuals with no data remaining before unnesting the track
### to avoid returning a null object. 

### Remove individuals with 0 rows first:
gps_unnest <- data_track_nest4 %>% 
  mutate(nrow = map_dbl(data, nrow)) %>% 
  filter(nrow > 0) %>% 
  # Unnest the data frame
  unnest(cols = data) %>% 
  # If you're worried about the unused factor levels, you
  # can redo the levels with 'factor()
  mutate(animal_id = factor(animal_id))

### Check that the levels are correct:
length(levels(gps_unnest$animal_id)) # should be as long as there are animals in your dataset

### ii. drop duplicates, fast steps, and fast roundtrips
gps_flagdrop <- gps_unnest %>% 
  filter(duplicate_ == FALSE & 
           fast_step_ == FALSE &
           fast_roundtrip_ == FALSE
  )




# VIII. Convert the dataset to a df and save as csv ----
gps_cleaned_df <- as.data.frame(gps_flagdrop)
names(gps_cleaned_df)[which(names(gps_cleaned_df) == "x_")] <- "longitude"
names(gps_cleaned_df)[which(names(gps_cleaned_df) == "y_")] <- "latitude"
write.csv(gps_cleaned_df, "ExampleScripts/Data/UDWR_Elk_BookCliffs_Cleaned.csv")




# IX. Build basic maps of seropos's and seroneg's ----
gps_sf <- st_as_sf(gps_cleaned_df,
                   coords = c("longitude", "latitude"),
                   crs = 4326)

# transform lsfs_sf to latlong
#gps_sf_latlong <- st_transform(gps_sf, crs = "+proj=longlat +datum=WGS84")

# extract bounding box
gps_bbox <- st_bbox(gps_sf)

# extract center coordinates for the domain of these elk points
cen_long <- (gps_bbox$xmin + gps_bbox$xmax) / 2
cen_lat <- (gps_bbox$ymin + gps_bbox$ymax) / 2

# play around a bit to get reasonable tiles at an appropriate zoom
m <- leaflet() %>%
  setView(lng=as.numeric(cen_long), 
          lat=as.numeric(cen_lat), zoom = 12) 
m %>% addProviderTiles("Stamen.Terrain") # you might have to do some nasty
# registration with google to access the tiles, but I think it's not necessary
# for Stamen.Terrain. 

# get fixes for females
fems <- subset(gps_cleaned_df, is.na(latitude) == F &
                 sex == "F")
# convert female fixes to sf object
fems_st <- st_as_sf(fems, coords = c("longitude", "latitude"), 
                    crs = 4326)
# get fixes for seronegative animals
males <- subset(gps_cleaned_df, is.na(latitude) == F &
                  sex == "M")
# convert seronegative fixes to sf object
males_st <- st_as_sf(males, coords = c("longitude", "latitude"), 
                     crs = 4326)

# build leaflet map
m <- leaflet(fems) %>% # first specify where to center the map (setView) and what zoom to use
  setView(lng=as.numeric(cen_long), 
          lat=as.numeric(cen_lat), zoom = 8) %>%
  addProviderTiles("Stamen.Terrain") # add a background map tile
m # print the map out so you can see whether the zoom looks right-ish.
m <- m %>% addCircles(lng = ~ longitude, # add location points for seronegative animals
                      lat = ~ latitude,
                      data = males,
                      color = "black") # these will NOT DISPLAY until you 
# reprint the "m" map (done several lines down)
m <- m %>% addCircles(lng = ~ longitude,
                      lat = ~ latitude,
                      data = fems,
                      color = "red")
# Now print the whole map. 
# this could take a second to render after printing (for my
# 24,000 points, it takes ~6 seconds)
m # prints map with contents specified above.

