## Introduction to Animal Movement Ecology
## K. Manlove, August 2023
##
## Intro to habitat selection analysis
## 

## This code is very lightly adapted from
#  Appendix A -- HSF Examples
#  John Fieberg, Johannes Signer, Brian Smith, Tal Avgar

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
elev_rast <- terra::rast("ExampleScripts/Data/BookCliffs_DerivedRasters/elev.tif")
slope_rast <- terra::rast("ExampleScripts/Data/BookCliffs_DerivedRasters/slope.tif")
aspect_rast <- terra::rast("ExampleScripts/Data/BookCliffs_DerivedRasters/aspect.tif")
tri_rast <- terra::rast("ExampleScripts/Data/BookCliffs_DerivedRasters/tri.tif")

# specify new crs (for UTMs). NOTE: epsg for UTMs in Utah is 26912. 
newcrs <- "+init=epsg:26912"
# project raster to UTMs using the project() function in terra
elev_rast_utm <- terra::project(elev_rast, newcrs)
slope_rast_utm <- terra::project(slope_rast, newcrs)
tri_rast_utm <- terra::project(tri_rast, newcrs)
aspect_rast_utm <- terra::project(aspect_rast, newcrs)
# rename elev_rast_utm to just be elev_rast for simplicity
names(elev_rast_utm) <- "elevation"
names(slope_rast_utm) <- "slope"
names(tri_rast_utm) <- "tri"
names(aspect_rast_utm) <- "aspect"
#### b. Landfire (https://landfire.gov/getdata.php) ----
vegtype_rast_01 <- terra::rast(x = "ExampleScripts/Data/BookCliffs_Landfire/landfire_5mT1bAOUHUjzZKf0FTGb/LF2022_EVH_230_CONUS/LC22_EVH_230.tif")
vegtype_rast_02 <- terra::rast(x = "ExampleScripts/Data/BookCliffs_Landfire/landfire_5mT1bAOUHUjzZKf0FTGb/LF2022_EVT_230_CONUS/LC22_EVT_230.tif")
# specify that the projection on the landfire data is the same as on the elev_rast
vegtype_rast_utm <- terra::project(vegtype_rast_01, crs(elev_rast_utm))
vegtype_rast_utm_02 <- terra::project(vegtype_rast_02, crs(elev_rast_utm))
s <- terra::sprc(vegtype_rast_utm, vegtype_rast_utm_02)
vegtype_rast <- terra::merge(s)
names(vegtype_rast) <- "vegtype"

# ADD THIS LINE!!!!
vegtype_resamp <- terra::resample(vegtype_rast, elev_rast_utm, method='bilinear')
vegtype_rast <- vegtype_resamp
names(vegtype_rast) <- "vegtype"

#### c. align rasters ----
# # set extent to use (here, extracted from the veg layer)
# extent_to_use <- terra::ext(vegtype_rast)
# # crop all other covariates to that extent
# elev_crop <- terra::crop(elev_rast_utm, extent_to_use)
# slope_crop <- terra::crop(slope_rast, extent_to_use)
# aspect_crop <- terra::crop(aspect_rast, extent_to_use)
# tri_crop <- terra::crop(tri_rast, extent_to_use)
# # resample vegtype so it matches attributes of DEM rasters
# vegtype_resamp <- terra::resample(vegtype_rast_utm, elev_crop, method='bilinear')

### iii. Append covariate data onto track and generate random points for one animal ----
# pull off just focal animal from track
animal_in <- "EL19F0068"
focal_trk <- trk %>%
  filter(id == animal_in) %>%
  arrange(t_)

summarize_sampling_rate(focal_trk)

## C. Format and generate random steps with environmental covariates attached ----
### i. generate random points ----
hsa_dat_focal_trk <- focal_trk %>% 
  random_points(n = length(focal_trk$x_) * 50) 
# I generated 50 randoms per used here.
# by default, points are sampled within an MCP for this animal.

### ii. append covariates at all used and available points ----
hsa_with_covs <- hsa_dat_focal_trk %>%
  extract_covariates(slope_rast_utm) %>%
  extract_covariates(tri_rast_utm) %>% 
  extract_covariates(aspect_rast_utm) %>% 
  extract_covariates(vegtype_rast) %>% 
  extract_covariates(elev_rast_utm)

### iii. scale covariates ----
hsa_covs_scaled <- hsa_with_covs %>%
  mutate(elev_rast = scale(elevation), # note that these names need to match the 
         # covariate names in hsa_with_covs (see names(hsa_with_covs) to print those)
         slope_rast = scale(slope)[, 1], 
         tri_rast = scale(tri)[, 1], , 
         low_tree = vegtype %in% c("Tree Height = 1 meter",  # recode landfire to more meaningful groupings
                                   "Tree Height = 2 meters",
                                   "Tree Height = 3 meters"),
         high_tree = vegtype %in% c("Tree Height = 4 meters",
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
         short_shrub = vegtype %in% c("Shrub Height = 0.2 meter",
                                      "Shrub Height = 0.3 meter",
                                      "Shrub Height = 0.4 meter",
                                      "Shrub Height = 0.5 meter",
                                      "Shrub Height = 0.6 meter",
                                      "Shrub Height = 0.7 meter",
                                      "Shrub Height = 0.8 meter",
                                      "Shrub Height = 0.9 meter",
                                      "Shrub Height = 1 meter"),
         tall_shrub = vegtype %in% c("Shrub Height = 1.1 meters",
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
         short_herb = vegtype %in% c("Herb Height = 0.1 meter",
                                     "Herb Height = 0.2 meter"),
         tall_herb = vegtype %in% c("Herb Height = 0.3 meter",
                                    "Herb Height = 0.4 meter"),
         weight = ifelse(case_, 1, 1e3))

names(hsa_covs_scaled)

### iv. add weights to available points to improve model's convergence on IPPP ----
hsa_covs_scaled$w <- ifelse(hsa_covs_scaled$case_, 1, 5000)


# II. Fit the logistic regression model to used and available points ----
# (stored in "case_"s)
## note that the fit is just a standard GLM fit. 
HSF_individ1 <- glm(case_ ~ elevation + slope + tri + low_tree + 
                      high_tree + short_shrub + tall_shrub + short_herb + tall_herb, 
                    data = hsa_covs_scaled, weight = w,
                    family = binomial(link = "logit"))


HSF_individ1 <- glm(case_ ~ elevation + slope + tri 
                    #+ low_tree + 
                    #  high_tree + short_shrub + tall_shrub + short_herb + tall_herb
                    , 
                    data = hsa_covs_scaled, weight = w,
                    family = binomial(link = "logit"))

summary(HSF_individ1)
confint(HSF_individ1)

# III. Model checking ----
## A. Assess model's randomized quantile residuals using DHARMa ----
library(DHARMa)
### i. simulate model residuals ----
# (this takes about 1 min on my machine)
elk_resids <- simulateResiduals(fittedModel = HSF_individ1,
                                plot = F, n = 100)

# to access the residuals:
residuals(elk_resids)
hist(residuals(elk_resids))
# these should look roughly uniform if the model is appropriate

# to generate the diagnostic plots:
# (these also take a minute to generate)
plot(elk_resids)

