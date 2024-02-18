## Introduction to Animal Movement Ecology
## K. Manlove, August 2023
##
## Section 3: Basic GIS operations in R
## 
##

# I. Load required libraries and data ----
## A. Libraries ----
library(sf)
library(sp)
library(terra)

## B. Response data ----
bookcliffs <- read.csv("ExampleScripts/Data/UDWR_Elk_BookCliffs.csv", header = T)
head(bookcliffs)
nrow(bookcliffs)

# convert dat_sm to spatial points using st_as_sf in sf
dat_sf <- sf::st_as_sf(bookcliffs, coords = c("longitude", "latitude"))
# apply crs to the new spatial points object
# NOTE: 4326 is the EPSG code for WGS84 lag-long. 
sf::st_crs(dat_sf) = 4326

## C. Bring in DEM and landfire products ----
### i. Bring in DEM ----
#### a. Read in and combine elevation tiles ----
elev_rast_01 <- terra::rast("ExampleScripts/Data/BookCliffs_DEM/USGS_13_n40w110_20220510.tif")
elev_rast_02 <- terra::rast("ExampleScripts/Data/BookCliffs_DEM/USGS_13_n40w111_20220510.tif")
# merge these two rasters into a common raster
s <- terra::sprc(elev_rast_01, elev_rast_02)
elev_rast <- terra::merge(s)
# save the merged raster to a new tif
writeRaster(elev_rast, "ExampleScripts/Data/BookCliffs_DerivedRasters/elev.tif")

# inspect the new raster
plot(elev_rast)

#### b. specify new crs (for UTMs). 
# NOTE: epsg for UTMs in Utah is 26912. 
newcrs <- "+init=epsg:26912"
# project raster to UTMs using the project() function in terra
elev_rast_utm <- terra::project(elev_rast, newcrs)
# rename elev_rast_utm to just be elev_rast for simplicity
names(elev_rast_utm) <- "elev_rast"
plot(elev_rast_utm)


### ii. Bring in landfire data ----
# https://landfire.gov/getdata.php
vegtype_rast_01 <- terra::rast(x = "ExampleScripts/Data/BookCliffs_Landfire/landfire_5mT1bAOUHUjzZKf0FTGb/LF2022_EVH_230_CONUS/LC22_EVH_230.tif")
vegtype_rast_02 <- terra::rast(x = "ExampleScripts/Data/BookCliffs_Landfire/landfire_5mT1bAOUHUjzZKf0FTGb/LF2022_EVT_230_CONUS/LC22_EVT_230.tif")
# specify that the projection on the landfire data is the same as on the elev_rast
vegtype_rast_utm <- terra::project(vegtype_rast_01, crs(elev_rast_utm))

plot(vegtype_rast_utm)



# III. Conduct spatial analyses to construct desired covariates ----
## A. Extract slopes
slope_rast <- terra::terrain(elev_rast_utm, v = "slope",
                             unit = "degrees", 
                             neighbors = 8)
writeRaster(slope_rast, "ExampleScripts/Data/BookCliffs_DerivedRasters/slope.tif")
## B. Extract aspects
aspect_rast <- terra::terrain(elev_rast_utm, v = "aspect",
                              unit = "degrees", 
                              neighbors = 8)
writeRaster(aspect_rast, "ExampleScripts/Data/BookCliffs_DerivedRasters/aspect.tif")

## C. Extract TRIs
tri_rast <- terra::terrain(elev_rast_utm, v = "TRI",
                           unit = "degrees", 
                           neighbors = 8)
writeRaster(tri_rast, "ExampleScripts/Data/BookCliffs_DerivedRasters/tri.tif")

# IV. Assemble covariate layers ----
## A. clip all DEM-derived rasters to extent sim to vegtype ----
extent_to_use <- terra::ext(vegtype_rast_utm)

elev_crop <- terra::crop(elev_rast_utm, extent_to_use)
slope_crop <- terra::crop(slope_rast, extent_to_use)
aspect_crop <- terra::crop(aspect_rast, extent_to_use)
tri_crop <- terra::crop(tri_rast, extent_to_use)

# clip all of them down to smaller rasters and save for future use. 
extent_to_use <- ext(600000, 670000, 4320000, 4420000)
elev_crop <- terra::crop(elev_rast_utm, extent_to_use)
slope_crop <- terra::crop(slope_rast_utm, extent_to_use)
aspect_crop <- terra::crop(aspect_rast_utm, extent_to_use)
tri_crop <- terra::crop(tri_rast_utm, extent_to_use)

writeRaster(elev_crop, "ExampleScripts/Data/BookCliffs_DerivedRasters/Smaller/elev_crop.tif")
writeRaster(slope_crop, "ExampleScripts/Data/BookCliffs_DerivedRasters/Smaller/slope_crop.tif")
writeRaster(aspect_crop, "ExampleScripts/Data/BookCliffs_DerivedRasters/Smaller/aspect_crop.tif")
writeRaster(tri_crop, "ExampleScripts/Data/BookCliffs_DerivedRasters/Smaller/tri_crop.tif")



### ii. resample vegtype so it matches attributes of DEM rasters ----
vegtype_resamp <- terra::resample(vegtype_rast_utm, elev_crop, method='bilinear')
extent_to_use <- ext(600000, 670000, 4320000, 4420000)
veg_crop <- terra::crop(vegtype_resamp, extent_to_use)
writeRaster(veg_crop, "ExampleScripts/Data/BookCliffs_DerivedRasters/Smaller/veg_crop.tif")

veg_crop <- terra::rast("ExampleScripts/Data/BookCliffs_DerivedRasters/Smaller/veg_crop.tif")
str(veg_crop)


### iii. build raster stack ----
cov_stack <- c(elev_crop, 
               slope_crop,
               aspect_crop,
               tri_crop,
               vegtype_resamp)


## E. Plot out all rasters ----
mat = matrix(c(1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3,
               0, 4, 4, 4, 4, 0, 0, 5, 5, 5, 5, 0), 
             nrow = 2, byrow = T)
pdf("Environmental_covariate_rasters.pdf", height = 8, width = 11)
layout(mat)
image(elev_crop, las = 1, xlab = "", ylab = "",
      main = "Elevation")
image(slope_crop, las = 1, col = terrain.colors(n = 100), xlab = "", ylab = "",
      main = "Slope")
image(aspect_crop, las = 1, xlab = "", ylab = "",
      main = "Aspect", col = topo.colors(n = 100))
image(tri_crop, las = 1, xlab = "", ylab = "",
      main = "Terrain ruggedness", col = topo.colors(n = 100))
image(vegtype_resamp, las = 1, xlab = "", ylab = "",
      main = "Vegetation type", col = cm.colors(n = 10))
dev.off()
