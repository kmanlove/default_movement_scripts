## Introduction to Animal Movement Ecology
## K. Manlove, August 2023
##
## Intro to step selection functions
## 

# 0. Overview of "competing rates" framing of multinomials ----
x_in <- seq(0, 15, by = 0.01)
exp_20.0 <- dexp(rate = 1, x = x_in)
exp_2.0 <- dexp(rate = .5, x = x_in)
exp_0.2 <- dexp(rate = .2, x = x_in)


# draw 100 batches of waiting times to event:
exp_20_draw <- rexp(n = 10000,  # Draw 100 "waiting times"
                    rate = 10) # "at rate = 20"
exp_2.0_draw <- rexp(n = 10000,  # Draw 100 "waiting times"
                     rate = .5) # "at rate = 20"
exp_0.2_draw <- rexp(n = 10000,  # Draw 100 "waiting times"
                     rate = .02) # "at rate = 20"

# cbind these together.
risk_sets <- cbind(exp_20_draw,
                   exp_2.0_draw,
                   exp_0.2_draw)
head(risk_sets)

# whichever event happened first in each risk set is the "event" we actually observed
observed_event <- apply(risk_sets, # apply to risk_sts
                        MARGIN = 1, # do same operation in every row
                        which.min) # here's the operation-to-do

# calculate proportion of occurrences of each event type
prob_of_event <- prop.table(table(observed_event))
prob_of_event


# plot out first set for illustration
par(mfrow = c(1, 3), mar = c(4, 6, 3, 1))
plot(exp_2.0 ~ x_in, type = "l", lwd = 3, 
     xlab = "Waiting time to event", 
     ylab = "Probability density", 
     las = 1, main = "Waiting time distributions")
lines(exp_0.2 ~ x_in, col = "blue", lwd = 3)
lines(exp_20.0 ~ x_in, col = "red", lwd = 3)
leg_text <- c("Event type 1", "Event type 2", "Event type 3")
legend("topright", leg_text, col = c("black", "blue", "red"),
       lwd = c(3, 3, 3), bty = "n")

plot(x = 0, y = 0, cex = 0, xlim = c(0, 3), ylim = c(0, 4),
     xlab = "Waiting time to event", ylab = "", las = 1,
     yaxt = "n", main = "Within one risk set")
segments(x0 = 0, x1 = risk_sets[1, 1], y0 = 1, y1 = 1, lwd = 3)
segments(x0 = 0, x1 = risk_sets[1, 2], y0 = 2, y1 = 2, lwd = 3, col = "blue")
segments(x0 = 0, x1 = risk_sets[1, 3], y0 = 3, y1 = 3, lwd = 3, col = "red")
abline(v = risk_sets[1, 1], lty = 2, lwd = 1)
axis(side = 2, at = c(1, 2, 3), labels = c("Event type 1", "Event type 2", "Event type 3"), las = 1)

#par(mfrow = c(1, 1), mar = c(4, 2, 3, 1))
plot(x = 0, y = 0, cex = 0, xlim = c(0, 3), ylim = c(0, 21),
     xlab = "Waiting time to event", ylab = "", las = 1,
     yaxt = "n", main = "Across many risk sets")
for(i in 1:20){
  segments(x0 = 0, x1 = risk_sets[i, 1], y0 = i + .3, y1 = i + .3, lwd = 2)
  segments(x0 = 0, x1 = risk_sets[i, 2], y0 = i + .5, y1 = i + .5, lwd = 2, col = "blue")
  segments(x0 = 0, x1 = risk_sets[i, 3], y0 = i + .8, y1 = i + .8, lwd = 2, col = "red")
  abline(v = risk_sets[i, which.min(risk_sets[i, ])], lty = 1, lwd = .5,
         col = c("black", "blue", "red")[which.min(risk_sets[i, ])])
}


# whichever event happened first in each risk set is the "event" we actually observed
observed_event <- apply(risk_sets, # apply to risk_sts
                        MARGIN = 1, # do same operation in every row
                        which.min) # here's the operation-to-do

# calculate proportion of occurrences of each event type
prob_of_event <- prop.table(table(observed_event))
prob_of_event



## This code is very lightly adapted from
#  Appendix B & C -- iSSA Examples
#  John Fieberg, Johannes Signer, Brian Smith, Tal Avgar
# https://conservancy.umn.edu/bitstream/handle/11299/218272/AppB_SSF_examples.html?sequence=26&isAllowed=y


# I. load required libraries and data ----
## A. Required libraries ----
library(tidyverse)
library(amt)
library(broom)
library(terra) # terra to manipulate covariate rasters

## B. Data ----
### i. Animal movement data ----
#### a. bring in data ----
bookcliffs <- read.csv("ExampleScripts/Data/UDWR_Elk_BookCliffs_Cleaned.csv", header = T)
bookcliffs$Range <- bookcliffs$cap_area
bookcliffs$uid <- bookcliffs$animal_id
bookcliffs$uid <- factor(bookcliffs$uid)
bookcliffs$dt <- lubridate::ymd_hms(as.character(bookcliffs$t_))
bookcliffs$season <- ifelse(month(bookcliffs$dt) %in% c(12, 1, 2, 3), "winter",
                            ifelse(month(bookcliffs$dt) %in% c(4, 5, 6), "spring",
                                   "summer-fall"))
table(bookcliffs$season)
bookcliffs_utm <- 32612  # epsg:32612 is for UTM Zone 12N
bookcliffs_latlong <- 4326

## check how many animals are present in the dataset
bookcliffs_animals <- length(levels(factor(bookcliffs$uid)))
bookcliffs_animals # should be 132

# check total relocations in the dataset
nrow(bookcliffs) # should be 700,001

# cut down to just the data from the first 20 animals
book_small <- subset(bookcliffs, uid %in% levels(factor(bookcliffs$uid))[1:20])
nrow(book_small) # 129,501 points for this set of animals
# reproject points to UTMs and append easting and northing onto the df.
# convert to sf object
book_small_sf <- sf::st_as_sf(book_small, coords = c("longitude", "latitude"))
# set coordinate reference system
sf::st_crs(book_small_sf) = 4326
# transform to UTMs
book_small_utmtransform <- st_transform(book_small_sf, crs = 32612)
# extract the coordinates from the transformed object into a dataframe
coords <- st_coordinates(book_small_utmtransform)
# append the easting and northing coords onto the original book_small dataframe
book_small$easting <- coords[, 1]
book_small$northing <- coords[, 2]
# get extent of book_small in utms
extent(book_small_utmtransform)

#### b. make amt track ----
trk <- book_small |> amt::make_track(.x = easting,
                                     .y = northing,
                                     .t = dt,
                                     id = uid,
                                     crs = bookcliffs_utm)

#### c. remove original dataframe to make room ----
rm(bookcliffs) # save space by removing the big df



### ii. Habitat rasters ----
#### a. DEMs ----
elev_rast <- terra::rast("ExampleScripts/Data/BookCliffs_DerivedRasters/Smaller/elev_crop.tif")
slope_rast <- terra::rast("ExampleScripts/Data/BookCliffs_DerivedRasters/Smaller/slope_crop.tif")
aspect_rast <- terra::rast("ExampleScripts/Data/BookCliffs_DerivedRasters/Smaller/aspect_crop.tif")
tri_rast <- terra::rast("ExampleScripts/Data/BookCliffs_DerivedRasters/Smaller/tri_crop.tif")

# rename elev_rast_utm to just be elev_rast for simplicity
names(elev_rast) <- "elevation"
names(slope_rast) <- "slope"
names(tri_rast) <- "tri"
names(aspect_rast) <- "aspect"

#### b. Landfire (https://landfire.gov/getdata.php) ----
veg_rast <- terra::rast("ExampleScripts/Data/BookCliffs_DerivedRasters/Smaller/veg_crop.tif")
levels(veg_rast)
names(veg_rast) <- "vegtype"

### iii. Append covariate data onto track and generate random points for one animal ----
# pull off just focal animal from track
animal_in <- "EL19F0068"
focal_trk <- trk %>%
  filter(id == animal_in) %>%
  arrange(t_)

summarize_sampling_rate(focal_trk)

## B. Format and generate random steps with environmental covariates attached ----
### i. resample the track to just (approximately) uniform timesteps ----
trk_common_ts <- focal_trk %>%
  track_resample(rate = hours(4), tolerance = minutes(20))

### ii. convert the track from points to steps ----
steps_common_ts <- trk_common_ts %>%
  steps_by_burst()

### iii. generate random steps ----
steps_with_randoms <- steps_common_ts %>% 
  random_steps(n_control = 20) # note that the random_steps argument is "n_control" instead of "n"
# I generated 20 random steps per realized step here. 
# by default, points are sampled within an MCP for this animal.

### iv. as an aside, take a look at the step length and turning angle distributions
# that random_steps worked from
sl_distr(steps_with_randoms) # details of the step length distribution
ta_distr(steps_with_randoms) # details of the turning angle distribution

### v. append covariates at all used and available points ----
ssf_dat_with_covs <- steps_with_randoms %>%
  extract_covariates(slope_rast, where = "both") %>% # the where = "both" argument tells
  # R to extract covariate values at both the beginning and the end of the step
  extract_covariates(tri_rast, where = "both") %>% 
  extract_covariates(aspect_rast, where = "both") %>% 
  extract_covariates(elev_rast, where = "both") %>%
  extract_covariates(veg_rast, where = "both")

# check the names on ssf_dat_with_covs -- those are the names you'll want to call below. 
names(ssf_dat_with_covs)

### vi. scale covariates ----
ssf_covs_scaled <- ssf_dat_with_covs %>%
  mutate(elev_start = scale(elevation_start), # note that these names need to match the 
         # covariate names in hsa_with_covs (see names(hsa_with_covs) to print those)
         elev_end = scale(elevation_end), 
         tri_start = scale(tri_start), 
         tri_end = scale(tri_end), 
         slope_start = scale(slope_start), 
         slope_end = scale(slope_end), 
         aspect_start = scale(aspect_start), 
         aspect_end = scale(aspect_end),
         low_tree_end = vegtype_end %in% c("Tree Height = 1 meter",  # recode landfire to more meaningful groupings
                                           "Tree Height = 2 meters",
                                           "Tree Height = 3 meters"),
         high_tree_end = vegtype_end %in% c("Tree Height = 4 meters",
                                            "Tree Height = 5 meters",
                                            "Tree Height = 6 meters",
                                            "Tree Height = 7 meters",
                                            "Tree Height = 8 meters",
                                            "Tree Height = 9 meters",
                                            "Tree Height = 10 meters",
                                            "Tree Height = 11 meters",
                                            "Tree Height = 12 meters",
                                            "Tree Height = 13 meters",
                                            "Tree Height = 14 meters",
                                            "Tree Height = 15 meters",
                                            "Tree Height = 16 meters",
                                            "Tree Height = 17 meters",
                                            "Tree Height = 18 meters",
                                            "Tree Height = 19 meters",
                                            "Tree Height = 20 meters",
                                            "Tree Height = 21meters"),
         short_shrub_end = vegtype_end %in% c("Shrub Height = 0.2 meter",
                                              "Shrub Height = 0.3 meter",
                                              "Shrub Height = 0.4 meter",
                                              "Shrub Height = 0.5 meter",
                                              "Shrub Height = 0.6 meter",
                                              "Shrub Height = 0.7 meter",
                                              "Shrub Height = 0.8 meter",
                                              "Shrub Height = 0.9 meter",
                                              "Shrub Height = 1 meter"),
         tall_shrub_end = vegtype_end %in% c("Shrub Height = 1.1 meters",
                                             "Shrub Height = 1.2 meters",
                                             "Shrub Height = 1.3 meters",
                                             "Shrub Height = 1.4 meters",
                                             "Shrub Height = 1.5 meters",
                                             "Shrub Height = 1.6 meters",
                                             "Shrub Height = 1.7 meters",
                                             "Shrub Height = 1.8 meters",
                                             "Shrub Height = 1.9 meters",
                                             "Shrub Height = 2.0 meters",
                                             "Shrub Height = 2.1 meters",
                                             "Shrub Height = 2.2 meters",
                                             "Shrub Height = 2.3 meters",
                                             "Shrub Height = 2.4 meters",
                                             "Shrub Height = 2.5 meters",
                                             "Shrub Height = 2.6 meters",
                                             "Shrub Height = 2.7 meters",
                                             "Shrub Height = 2.8 meters",
                                             "Shrub Height = 2.9 meters",
                                             "Shrub Height >= 3.0 meters"),
         short_herb_end = vegtype_end %in% c("Herb Height = 0.1 meter",
                                             "Herb Height = 0.2 meter"),
         tall_herb_end = vegtype_end %in% c("Herb Height = 0.3 meter",
                                            "Herb Height = 0.4 meter"),
         
         cos_ta_ = cos(ta_), # also need to add on cos of the turning angle
         # (which amt built for us when we built steps_common_ts)
         # and log of the step length
         log_sl_ = log(sl_))

# do a little inspection, especially on the categorical predictors.
names(ssf_covs_scaled)
ssf_covs_scaled$elev_start[1:100]
ssf_covs_scaled$high_tree_end[1:10]

### vii. remove any steps with missing turn angles ----
ssf_covs_final <- ssf_covs_scaled %>%
  filter(!is.na(ta_))

names(ssf_covs_final)

# II. Fit an iSSF model with environmentally-determined habitat selection  ----
# but environmentally-independent movements
issf_hab_seln_with_constant_movement <- ssf_covs_final %>% 
  fit_issf(case_ ~ slope_end + elev_end + tri_end + # these covs drive the end (= hab seln) process
             sl_ + log_sl_ + cos_ta_ + # need to include step length, log(step length)
             # and cos(turning angle) alongside covariates in the model.
             strata(step_id_), # this tells fit_issf to fit a conditional logit
           # to the data arising from each step. 
           model = TRUE)

summary(issf_hab_seln_with_constant_movement)


# III. Fit an iSSF model with environmentally-determined habitat selection & movement ----
issf_hab_seln_environ_move <- ssf_covs_final %>% 
  fit_issf(case_ ~ slope_end + elev_end + tri_end + # these covs drive the end (= hab seln) process
             slope_start:(sl_ + log_sl_ + cos_ta_) + # for environment to effect movement, need
             # environmental attributes to interact with movement attributes (sl, log_sl, cos_ta)
             sl_ + log_sl_ + cos_ta_ + # need to include step length, log(step length)
             # and cos(turning angle) alongside covariates in the model.
             strata(step_id_), # this tells fit_issf to fit a conditional logit
           # to the data arising from each step. 
           model = TRUE)

summary(issf_hab_seln_environ_move)



