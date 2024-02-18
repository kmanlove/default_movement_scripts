## Introduction to Animal Movement Ecology
## K. Manlove, August 2023
## Lauren Ricci (2018-2019)
## Grete Wilson-Henjum (2021-2022)
## Josh O'Brien (2019-2022)
##
## Building range and occurrence distributions in R
## 

# I. Libraries and data set-up ----
## A. Libraries ----
library(terra)
library(adehabitatHR)
library(sf)
library(lubridate)
library(ctmm)
library(raster)

## B. Data set-up ----

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

## Your dataframe should contain the following fields:
### Range  * a variable containing a herd name. If your data are just from one herd, you
###        * can ignore this field and read your original data in in lieu of the "zion" df.
### uid OR ID *contains animal IDs (adjust code accordingly depending on how this is named)
### a date-time field that can be appropriately extracted using lubridate
### x, y   * in two separately columns (these could be UTMs or Lat-long; 
###        * this script is set up for UTMs from zone 12N labeled EAST and NORTH, which I relabel below.
###        * You'll need to change the zion_crs field to reflect your coordinate reference system
###        * if this is NOT your zone, or if your data are in lat-long. 
### season * a categorical field indicating what season each point falls into (this
###        * could be built on the fly using your date field)

# ## read in the sample gps dataset (it's big)
# #allranges <- read.csv("SampleData/gps_allranges_new_20220408.csv", header = T)
# ## subset down to just the Zion data
# #zion <- subset(allranges, Range == "zion")
# ## check the dataset's size to be sure you've got it all:
# nrow(zion) # should be 154,630
# ## specify a coordinate reference system for these data
# zion_crs <- 32612  # epsg:32612 is for UTM Zone 12N
# ## check how many animals are present in the dataset
# zion_animals <- length(levels(factor(zion$uid)))
# zion_animals # should be 22
# ## Set up the date-time fields 
# ###        * The Zion data have a PACDATE field and a separate TIME field, which
# ###        * I'm pasting together and passing to lubridate:: mdy_hms to make dt.
# ###        * You may need to use a different set of functions
# ###        * here depending on how your dates and times are formatted.
# zion$dt <- lubridate::mdy_hms(paste(zion$PACDATE, " ", zion$TIME, sep = ""))
# 
# ## rename the ID field "uid" 
# ## (this is the field name that some of the packages look for for the animal's ID)
# names(zion)[which(names(zion) == "ID")] <- "uid"
# ## convert uid field to be a factor
# zion$uid <- factor(zion$uid)
# 
# # add "x" and "y" columns that are copies of EAST and NORTH
# zion$x <- as.numeric(as.character(zion$EAST))
# zion$y <- as.numeric(as.character(zion$NORTH))
# 
# zion <- subset(zion, x <= 360000 & y <- 4150000)

## C. Set up data and execute functions for each animal in the dataset ----
### i. data set-up ----
test_unique_ids <- levels(factor(book_small$uid))
individ_list <- vector("list", length(test_unique_ids))
crs_in <- bookcliffs_utm

# Extract bounding box on raster for whole herd (building all individual HRs on the same raster
# ensures that we'll easily be able to stack the individual rasters up to build a herd aggregate)
# I add 5K of extra space around edges for kde and other hr estimators with Gaussian smoothing.
## 1. convert dataframe to sf object
data_in_sf <- sf::st_as_sf(book_small, coords = c("longitude", "latitude"), crs = bookcliffs_latlong)
data_in_utms <- sf::st_transform(data_in_sf, crs = bookcliffs_utm) # reproject data to UTMs
data_in_sp <- as_Spatial(data_in_utms) # convert to SpatialPointsDataFrame

geom <- st_coordinates(data_in_utms)
book_small$x <- geom[, 1]
book_small$y <- geom[, 2]


## 2. build a raster defined by the bounds of that sf object
### all the homeranges we build will sit on top of this raster. This keeps them all
### having the same projection, extent, etc. 
grand_raster <- raster::raster(xmn = xmin(sf::as_Spatial(data_in_utms))-10000, 
                               # (gets minimum x coordinate, and subtracts another 5000m to set the box's southern bound)
                               xmx = xmax(sf::as_Spatial(data_in_utms))+10000, 
                               # (gets maximum x coordinate, and adds another 5000m to set the box's northern bound)
                               ymn = ymin(sf::as_Spatial(data_in_utms))-10000, 
                               ymx = ymax(sf::as_Spatial(data_in_utms))+10000, 
                               res = c(200, 200), # sets the resolution of the raster to 100m x 100m... 
                               # runs faster if resolution is lower. 
                               crs = bookcliffs_utm
)
# Define coordinate reference system for grand_raster to be crs for this herd. 
crs(grand_raster) <- paste("+proj=utm +init=epsg:", bookcliffs_utm, sep = "")

# Convert grand raster to spatial pixels, spatial grid, and spdf. 
# We do this because different functions require this object to be formatted in different ways. 
grand_pixels <- as(grand_raster, "SpatialPixelsDataFrame") 
grand_grid <- as(grand_raster, "SpatialGridDataFrame")          # throws warning; don't worry. 
grand_pixels_spdf <- as(grand_pixels, "SpatialPixelsDataFrame") # throws warning; don't worry.



# II. Code to build a single home range of each type for one animal ----
## A. set up data for one animal ----
# cut down to just the first animal in book_small
individ_in <- subset(book_small, uid == levels(factor(book_small$uid))[1], 
                     select = c("x", "y", "uid"))
# convert that animal's data to an sf object
individ_in <- sf::st_as_sf(individ_in, coords = c("x", "y"), crs = bookcliffs_utm)
# convert from sf to SpatialPoints for adehabitatHR
individ_in <- sf::as_Spatial(individ_in)
# reset individual ID to be a factor (as character, it crashes kernelUD below)
individ_in$uid <- factor(individ_in$uid)

## B. Fit MCP ----
mcp_fit_95 <- mcp(xy = individ_in, 
                  percent = 95, 
                  unin = c("m", "km"),
                  unout = c("ha", "km2", "m2"))
mcp_fit_50 <- mcp(xy = individ_in, 
                  percent=50, 
                  unin = c("m", "km"),
                  unout = c("ha", "km2", "m2"))
# extract the area within the 95% MCP
mcp_95Perc_area <- as.numeric(mcp.area(individ_in, percent = 95,
                                       unout = c("ha"), 
                                       plotit = F)[1])
mcp_95Perc_area
# extract area within the 50% MCP
mcp_50Perc_area <- as.numeric(mcp.area(individ_in, percent = 50,
                                       unout = c("ha"), plotit = F)[1])
mcp_95Perc_area


## C. Fit KDE ----
# fit a conventional kde to the points
kde_fit <- adehabitatHR::kernelUD(xy = individ_in,
                                  h = "href", # use the "href" bandwdith smoothing procedure
                                  kern = "epa", # use an Epanechnikov kernel
                                  grid = grand_pixels) # map the KDE out on grand_pixels (which is set up to be applicable to all animals)

# extract the 95% isopleth -- DO NOT WORRY ABOUT CRS WARNINGS.
kde_95Perc_area <- as.numeric(kernel.area(kde_fit,  # pass in KDE fit
                                          percent = 95, # set isopleth of interest
                                          unout = "ha")) # set units for area (here, hectares)
kde_95Perc_area
# extract the 50% isopleth
kde_50Perc_area <- as.numeric(kernel.area(kde_fit, 
                                          percent = 50, unout = "ha"))
kde_50Perc_area
# get boundaries for 95% and 50% isopleth
kde_95Perc_boundary <- getverticeshr(kde_fit, percent = 95)
kde_50Perc_boundary <- getverticeshr(kde_fit, percent = 50)

plot(kde_95Perc_boundary)
lines(kde_50Perc_boundary, col = "grey60")

## D. Fit Brownian bridge ----
# reset data to just be a df of this individual's points
individ_in <- subset(book_small, uid == levels(factor(book_small$uid))[1])

# build ltraj object (NOTE: I'm using UTMs here b/c it's easier to set the spatial
# imprecision in meters than degrees)
individ_traj <- adehabitatLT::as.ltraj(xy = individ_in[, c("x", "y")], 
                                       date = individ_in$dt, 
                                       id = factor(individ_in$uid), 
                                       proj4string = CRS(paste("+proj=utm +init=epsg:", bookcliffs_utm, sep = "")))
# plot the ltraj to be sure it looks alright
plot(individ_traj)

## fit Brownian Bridge
sig2 <- 10 # specify approximate spatial imprecision of points.
# estimate animal's speed parameter (sig1) using sig2 and the liker function 
sig1_est <- adehabitatHR::liker(individ_traj, 
                                sig2 = sig2, 
                                rangesig1 = c(0, 100), 
                                plotit = F)

sig1_val <- lapply(sig1_est, getElement, "sig1") %>% 
  unlist() %>% 
  unname()
sig1_val

# NOTE: THIS NEXT LINE TAKES APPROXIMATELY 3 MINS ON MY MACHINE
bbmm <- adehabitatHR::kernelbb(individ_traj, 
                               sig1 = sig1_val, 
                               sig2 = sig2
                               #,
                               #grid = grand_pixels
)

summary(bbmm)
bb_95 <- adehabitatHR::getverticeshr(bbmm, percent = 95)
bb_area_95 <- raster::area(bb_95) / 10000

plot(bbmm) # note, this is a super-coarse grid; could do better with tinkering. 


## E. Fit LoCoH ----
individ_in <- subset(book_small, uid == levels(factor(book_small$uid))[1], 
                     select = c("x", "y", "uid"))
# convert that animal's data to an sf object
individ_in <- sf::st_as_sf(individ_in, coords = c("x", "y"), crs = bookcliffs_utm)
# convert from sf to SpatialPoints for adehabitatHR
individ_in <- sf::as_Spatial(individ_in)
# reset individual ID to be a factor (as character, it crashes kernelUD below)
individ_in$uid <- factor(individ_in$uid)
k_in <- sqrt(nrow(individ_in)) # this is the value of k_in suggested in Getz. However,
# sqrt(# rows) is really big in this case (~53) and is messing with the model's ability to fit appropriately.
# So, for now, I'm not actually using k_in, just overriding and using k = 5, which is the default. 
## NEXT LINE IS SLOW.
locoh_fit <- adehabitatHR::LoCoH.k(individ_in, 
                                   k = 5, # for now, I'm using k instead of a. 
                                   # lower values of "a" lead to fewer orphaned hole errors. 
                                   duplicates = "random")
crs(locoh_fit) <- paste("+proj=utm +init=epsg:", bookcliffs_utm, sep = "")
locoh_area_95 <- adehabitatHR::LoCoH.k.area(individ_in,
                                            krange = k_in,
                                            unin = "m",
                                            unout = "km2")
plot(locoh_fit)


## F. Fit AKDE ----
### i. reset the data so it works with ctmm ----
# start back at just the first animal's locations from book_small as a df
individ_in_df <- subset(book_small, uid == levels(factor(book_small$uid))[1])
nrow(individ_in_df)
head(individ_in_df)

# reformat the data to match movebank naming conventions
# to convert to telemetry object in ctmm.
# names should be:
#   individual.local.identifier
#   timestamp
#   location.long
#   location.lat
individ_in_df$individual.local.identifier <- individ_in_df$animal_id
individ_in_df$timestamp <- individ_in_df$datePOSIX
individ_in_df$location.long <- individ_in_df$Longitude
individ_in_df$location.lat <- individ_in_df$Latitude

# animal_telem <- ctmm::as.telemetry(individ_in_df,
#                                    drop = FALSE) 

# build ou model to extract estimate of tau (time to independence)
m_ou <- ctmm(omega = 23, tau = 50) 
# specify model starting values
# tau = 50d is our best starting bet at time to statistical independence (TTSI)
### (this is similar to home range crossing time) based on past experience. 
## omega = 23km^2 is completely made up.
# define individ_telem for just one animal. I'm cutting dop, since ctmm will use it in its construction
# of the telemetry object if it's included, and I've already cleaned these data. 
individ <- subset(individ_in_df, animal_id == levels(factor(individ_in_df$animal_id))[1])
individ <- individ[, !(names(individ) == "dop")]
# convert to telemetry object
individ_telem <- ctmm::as.telemetry(individ,
                                    drop = FALSE) 
#individ_ctmm_fit_ou <- ctmm::ctmm.fit(individ_telem[[1]], m_ou) # fit the ou model
# extract relevant bits of output. 
# ou_tau_est <- summary(individ_ctmm_fit_ou)$CI[2, 2]
# ou_tau_lb <- summary(individ_ctmm_fit_ou)$CI[2, 1]
# ou_tau_ub <- summary(individ_ctmm_fit_ou)$CI[2, 3]
# ou_area_est <- summary(individ_ctmm_fit_ou)$CI[1, 2]
# ou_area_lb <- summary(individ_ctmm_fit_ou)$CI[1, 1]
# ou_area_ub <- summary(individ_ctmm_fit_ou)$CI[1, 3]

# now ID and fit the best autocorrelated model.
# see: https://cran.r-project.org/web/packages/ctmm/vignettes/akde.html
# establish OU starting conditions for this individual
individ_ctmm_ouf <- ctmm::ctmm.guess(individ_telem[[1]], interactive = FALSE) # automated model guess
# fit OU CTMM for this individual. 
# NOTE: THE NEXT LINE IS SLOWISH (~ 90 seconds on my machine)
individ_ctmm_fit <- ctmm::ctmm.fit(individ_telem[[1]], individ_ctmm_ouf) # automated model guess
# AKDE fit with automated guess
individ_akde_fit <- ctmm::akde(individ_telem[[1]], individ_ctmm_fit)
plot(individ_akde_fit)
path_in <- "ExampleScripts/Output/rangeoccurrence/akde/"
individ_akde_raster <- ctmm::writeRaster(individ_akde_fit, 
                                         filename = paste(path_in, "akde_fit_raster.grd", sep = ""),
                                         DF = "PDF")
saveRDS(individ_akde_fit, file = paste(path_in, "akde_fit.rds", sep = ""))
individ_summary <- summary(individ_akde_fit)
individ_summary
akde_area_95 <- individ_summary$CI[2] * 100
akde_95 <- SpatialPolygonsDataFrame.UD(individ_akde_fit, level.UD = 0.95)



# III. Occurrence overlaps in amt ----
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8048401/
library(amt)
book_small$uid_num <- as.numeric(as.factor(book_small$uid))

book_small_dat <- book_small %>% 
  # mutate(sex = substr(id, 1, 1)) %>% 
  make_track(longitude, latitude, dt, 
             id = uid_num,
             crs = bookcliffs_latlong)

book_small_trast <- make_trast(book_small_dat, 
                               res = 0.00000001,
                               factor = 1.1) 
# I had to increase the value of "factor" to get the KDE to build properly. 
# The documentation on what "factor" is doing is pretty scant, however. 

# nest the data
book_small_nest <- book_small_dat %>% 
  nest(data = c(x_, y_, t_)) 

# add hr columns
hr1 <- book_small_nest %>% 
  mutate( 
    hr_mcp = map(data, hr_mcp, trast = book_small_trast),
    hr_kde = map(data, hr_kde, trast = book_small_trast)
  )  


hr_overlap_mat <- matrix(NA, 
                         nrow = length(levels(factor(book_small$uid))),
                         ncol = length(levels(factor(book_small$uid))))
for(i in 1:(length(levels(factor(book_small$uid_num))) - 1)){
  for(j in (i + 1): length(levels(factor(book_small$uid_num)))){
    ## for MCPs
    # hr_overlap_ij <- hr_overlap(hr1$hr_mcp[[i]], hr1$hr_mcp[[j]], type = "hr")$overlap
    ## for KDEs
    hr_overlap_ij <- hr_overlap(hr1$hr_kde[[i]], 
                                hr1$hr_kde[[j]], 
                                type = "hr")$overlap
    hr_overlap_mat[i, j] <- hr_overlap_mat[j, i] <- hr_overlap_ij
  }
  print(i)
}

hr_overlap_mat
diag(hr_overlap_mat) <- rep(0, nrow(hr_overlap_mat)) # removes selfing loops in the graph

# bring in igraph
library(igraph)

# convert my matrix to a graph object in igraph
elk_graph <- graph_from_adjacency_matrix(hr_overlap_mat, 
                                         weighted = T) # specify graph has edge weights, not just binary edges

# specifying coordinates for a plot of the graph (using the Fruchterman-Reingold algorithm)
elk_coords_fr <- layout.fruchterman.reingold(elk_graph) # specify a graph layout that pulls

# nodes with stronger overlap toward one another in the graph plot
plot(elk_graph, # plots the graph
     layout = elk_coords_fr) # tells plot to use the fr coordinates from above

# run community detection on the graph
# Kezia's default mode is to use walktrap with 4 steps
elk_walktrap_clusters <- cluster_walktrap(elk_graph,
                                          weights = E(elk_graph)$weight,
                                          steps = 4)
# examine the resulting communities object
# membership slot tells us the community that each node belongs to
memberships <- elk_walktrap_clusters$membership
table(memberships)

# from here, extract all MCPs (or points, or whatever G-space object we have)
# for all animals within a cluster, and then do that for all clusters, and then map the
# clusters out in G-space

# extract members of one cluster
grp_1 <- which(memberships == 1)

# extract HRs for those animals
grp_1_tracks <- book_small_nest$data[which(book_small_nest$id %in% grp_1)]
class(grp_1_tracks)

# plot out HRs for this set in G-space
library(leaflet)
library(mapview)
gps_bbox <- amt::bbox(book_small_dat,
                      spatial = F) # if spatial != F, amt::bbox gives you back a polygon. 

# extract center coordinates for this dataset (here, book_small)
cen_long <- (gps_bbox$xmin + gps_bbox$xmax) / 2
cen_lat <- (gps_bbox$ymin + gps_bbox$ymax) / 2

# build leaflet map
m <- leaflet() %>% # first specify where to center the map (setView) and what zoom to use
  setView(lng=as.numeric(cen_long), 
          lat=as.numeric(cen_lat), zoom = 10) %>%
  addProviderTiles("Stamen.Terrain") # add a background map tile
# you might have to do some nasty registration with google to access the tiles, 
# but I think it's not necessary for Stamen.Terrain. 
m # print the map out so you can see whether the zoom looks right-ish.

# need to project the mcp back to lat-long. this is annoying. 

# version of maps that adds the home range polygons ----
m <- leaflet(hr1$hr_mcp[[1]]$mcp) %>% # first specify where to center the map (setView) and what zoom to use
  setView(lng=as.numeric(cen_long), 
          lat=as.numeric(cen_lat), zoom = 10) %>%
  addProviderTiles("Stamen.Terrain") %>%
  addPolygons(color = "#444444", weight = 1, smoothFactor = 0.5,
              opacity = 1.0, fillOpacity = 0.5,
              highlightOptions = highlightOptions(color = "white", weight = 2,
                                                  bringToFront = TRUE)) # add a background map tile
# you might have to do some nasty registration with google to access the tiles, 
# but I think it's not necessary for Stamen.Terrain. 
m # print the map out so you can see whether the zoom looks right-ish.

