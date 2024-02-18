## Introduction to Animal Movement Ecology
## K. Manlove, August 2023
##
## Building and analyzing ctmms
## 

# I. Data preparation and construction of telemetry object(s) ----
## A. load data and libraries ----
library(ctmm)
library(sp)
library(sf)
library(lubridate)

# read in the csv of points
bookcliffs <- read.csv("ExampleScripts/Data/UDWR_Elk_BookCliffs_Cleaned.csv", header = T)
# convert points to an sf object (similar to a "spatial points" object in a GIS)
dat_sf <- sf::st_as_sf(bookcliffs, coords = c("longitude", "latitude"))

# apply crs to the new spatial points object
# NOTE: 4326 is the EPSG code for WGS84 lag-long. 
sf::st_crs(dat_sf) = 4326

# extract (lat-long) coordinates and bind them onto the dataframe in named columns
coords_wgs84 <- sf::st_coordinates(dat_sf)
dat_sf$Latitude <- coords_wgs84[, 2]
dat_sf$Longitude <- coords_wgs84[, 1]

# ctmm would like the data in UTMs, not lat-long, so reproject. 
# NOTE: epsg for UTMs in Utah is 26912. 
newcrs <- "+init=epsg:26912"
# project raster to UTMs using the st_transform() function in sf
dat_utm <- sf::st_transform(dat_sf, newcrs)

# extract (UTM) coordinates and append them to dataframe in named columns
coords <- st_coordinates(dat_utm)
dat_utm$East <- coords[, 1]
dat_utm$North <- coords[, 2]

# convert the data to a POSIX object using lubridate
dat_utm$datePOSIX <- lubridate::ymd_hms(as.character(dat_utm$t_))

# cut down to just the cow elk
cows <- subset(dat_utm, sex == "F")
cows$animal_id <- factor(cows$animal_id)
cows_in <- as.character(levels(factor(cows$animal_id)))

## B. Construct & prep data for each animal ----
for(anim in 1:length(cows_in)){
  # subset full dataset down to just the focal cow
  animalIn <- cows_in[anim]
  # prep data for moveHMM
  animal_pts <- subset(cows, animal_id == animalIn)
  # strip off the geometry (the geometry gets written as two columns if you don't
  # drop it, and messes up the names of the columns when you read the csv back
  # in. Since we've already appended the geometry components as named columns, 
  # I'm just stripping them off here). 
  animal_pts_nogeom <- sf::st_drop_geometry(animal_pts)
  write.csv(animal_pts_nogeom, file = paste("ExampleScripts/Data/BookCliffs_Fixes/ByAnimal/", animal_pts$animal_id[1],  ".csv", sep = ""))
  print(anim)
}

# read in all animal-specific datasets and bind them together 
# (this will just regenerate the oridinal dat_sf)
animal_csvs <- list.files("ExampleScripts/Data/BookCliffs_Fixes/ByAnimal/")
animal_list <- vector("list", length(animal_csvs))
for(i in 1:length(animal_csvs)){
  animal_list[[i]] <- read.csv(paste("ExampleScripts/Data/BookCliffs_Fixes/ByAnimal/", animal_csvs[i], sep = ""))
}

animal_full <- do.call("rbind", animal_list)
dim(animal_full)
head(animal_full)

## C. reformat the data to match movebank naming conventions ----
# to convert to telemetry object in ctmm.
# names should be:
#   individual.local.identifier
#   timestamp
#   location.long
#   location.lat
names(animal_full)
animal_full$individual.local.identifier <- animal_full$animal_id
animal_full$timestamp <- animal_full$datePOSIX
animal_full$location.long <- animal_full$Longitude
animal_full$location.lat <- animal_full$Latitude

# pull out dop column, since ctmm::as.telemetry() tried to interpret it
# (and we've already dealt with bounce in the datacleaning protocol in script 02)

animal_full_names <- names(animal_full)
animal_full <- subset(animal_full, 
                      select = animal_full_names[-which(animal_full_names == "dop")])

# convert to telemetry object
animal_telem <- ctmm::as.telemetry(animal_full,
                                   drop = FALSE) 
# (drop = FALSE specifies that we're building a list of telemetry objects, one per animal)
class(animal_telem)
plot(animal_telem[[1]]) # plot the first animal's telemetry object

## D. check for outliers ----
outliers_1 <- ctmm::outlie(animal_telem[[1]], plot = TRUE, by = "d")
plot(outliers_1)


# II. Construction of visual diagnostics ----
## A. Empirical variograms ----
variograms <- lapply(animal_telem, variogram)

par(mfrow = c(4, 6), 
    oma = c(0, 0, 0, 0),# I'm adjusting margins here so I can
    # see the structure of the semivariograms better
    mar = c(3, 4, 0, 0)) # same for mar (inner margins)
for(i in 1:24){ # plot the variograms for the first 24 animals
  plot(variograms[[i]])
}

# zoom in on the very beginning of the variogram(s)
# to see whether there's also autocorrelation in
# velocity (in which case this curve would be sinusoidal in 
# shape). 
par(mfrow = c(1, 2))
plot(variograms[[1]])
# In this case, the shape looks pretty linear.
plot(variograms[[1]], fraction = 0.005)

# things to note on the semivariogram plot:
# 1) WHERE THE LINE CROSSES THE Y-AXIS (the "nugget"): 
#    note that the curves all start off at ~0 at
#    time = 0.  If there were consistent error in locations,. 
#    the line wouldn't drop all the way to 0 on the y-axis. 
#    That would be called a "nugget effect", but we see no
#    evidence of one here.  (as expected -- spatial error
#    is likely quite low in these collars)
# 2) WHERE THE LINE IS STEADILY INCREASING (region of autocorrelation):
#    The region of increase (before the curves flatten out)
#    is the region of (most severe) autocorrelation. For these
#    animals, that period of releant autocorrelation bounces
#    around quite a bit.  (for the 3rd and 4th animals, it might
#    be ~6 months, but for the 5th and 6th, it might be more like
#    1 month). 
# 3) IF/WHERE THE LINE ASYMPTOTES:
#    The line will hit an asymptote for range residents, and the
#    height of that line will be proportional to home range size. 
#    Here, most animals look like range residents, but there are a
#    few exceptions (for example, the 3rd and 23rd animals)

# We can empirically fit the variograms using the variogram.fit
# function (I _think_ the fits start at the mle, and the slides
# let you move things around from there. )
### i. Full variograms ----
variogram_fits <- lapply(variograms, variogram.fit)

### ii. Zoom in on just the early steps ----
lapply_vario_early <- function(x){
  variogram.fit(x, fraction = 0.005) # fraction = 0.005 restricts the plotted timeseries
  # to just the first half a percent of all points
}
variogram_fits_earlysteps <- lapply(variograms, lapply_vario_early)


# III. Fitting and comparing CTMMs for one animal  ----
#ctmm_fit_1 <- ctmm.fit(animal_telem[[1]],
#                       CTMM = ctmm())

## A. Identify initial conditions for ctmm for first animal ----
# these will be starting parameters for the ctmm fit a couple of lines further down. 
ctmm_guess_1 <- ctmm.guess(animal_telem[[1]], interactive = F)

## B. Bit the ctmm models for animal 1, applying initial conditions from ctmm_guess_1. ----
# (here, I'm also tracking runtime in my code.)
start <- proc.time() # proc.time() stores the time at which this line is executed.
# I'm using it to keep track of how long it takes to fit the ctmm.
ctmm_fit_1 <- ctmm.fit(animal_telem[[1]], 
                       CTMM = ctmm_guess_1,
                       COV = T)
end <- proc.time() # record the time at which this line executes
end - start # take the difference between end and start -- this is time elapsed
# for me, the fit for the first animal took about a minute. 

## C. Examine output ----
ctmm_fit_1 # many slots....
summary(ctmm_fit_1) # a little better organized
ctmm_fit_1$tau   # temporal autocorrelation
ctmm_fit_1$sigma # isotropy

## D. Rinse and repeat using ctmm.select() ----
# now, use ctmm.select to select among a set of candidate movement models for 
# one animal (start, end, and runtime are all helper calls that let me 
# track runtime in the code)
start <- proc.time()
ctmm_select_1 <- ctmm.select(animal_telem[[1]], 
                             CTMM = ctmm_guess_1, 
                             verbose = T)
end <- proc.time()
runtime <- end - start
runtime # 3.9 minutes to run.

ctmm_select_1
summary(ctmm_select_1)
plot(ctmm_select_1)

# look at the models in isolation
OUF_aniso_1 <- ctmm_select_1[[1]]
OUF_iso_1 <- ctmm_select_1[[2]]
OU_aniso_1 <- ctmm_select_1[[3]]
OUf_aniso_1 <- ctmm_select_1[[4]]

## E. Plot semivariograms ----
#pdf("MovementModels/ElkCtmm_variogramExamples.pdf", height = 7, width = 10)
par(mfrow = c(2, 4))
plot(variograms[[1]], main = "OUF (anisotropic)",
     CTMM = OUF_aniso_1, col.CTMM = rgb(0, 0, 1, alpha = .8))
plot(variograms[[1]], main = "OUF (isotropic)",
     CTMM = OUF_iso_1, col.CTMM = rgb(0, 0, 1, alpha = .8))
plot(variograms[[1]], main = "OU (anisotropic)",
     CTMM = OU_aniso_1, col.CTMM = rgb(0, 0, 1, alpha = .8))
plot(variograms[[1]], main = "OUf (anisotropic)",
     CTMM = OUf_aniso_1, col.CTMM = rgb(0, 0, 1, alpha = .8))

plot(variograms[[1]], main = "OUF (anisotropic)", CTMM = OUF_aniso_1, col.CTMM = rgb(0, 0, 1, alpha = .8),
     fraction = 0.005)
plot(variograms[[1]], main = "OUF (isotropic)", CTMM = OUF_iso_1, col.CTMM = rgb(0, 0, 1, alpha = .8),
     fraction = 0.005)
plot(variograms[[1]], main = "OU (anisotropic)", CTMM = OU_aniso_1, col.CTMM = rgb(0, 0, 1, alpha = .8),
     fraction = 0.005)
plot(variograms[[1]], main = "OUf (anisotropic)", CTMM = OUf_aniso_1, col.CTMM = rgb(0, 0, 1, alpha = .8),
     fraction = 0.005)

#dev.off()


# IV. Fitting and summarizing ctmms for all animals in dataset ----
## A. build and store ctmms for all animals ----
length(animal_telem) # 119 elk in this dataset. 
for(i in 1:length(animal_telem)){
  start <- proc.time() # proc.time() stores the time at which this line is executed.
  
  # fitting the ctmm for animal 23 (who doesn't have much of a sill)
  # identify initial conditions for ctmm for first animal
  ctmm_guess_i <- ctmm.guess(animal_telem[[i]], interactive = F)
  
  # fit the ctmm models for animal 1, applying
  # initial conditions from ctmm_guess_1.
  ctmm_fit_i <- ctmm.fit(animal_telem[[i]], 
                         CTMM = ctmm_guess_i,
                         COV = T)
  saveRDS(ctmm_fit_i, paste("ExampleScripts/Output/ctmms/ctmm_guess/", names(animal_telem)[i], ".rds", sep = ""))
  
  
  # select among ctmm models for one animal
  # (start, end, and runtime are all
  # helper calls that let me track runtime in my code)
  ctmm_select_i <- ctmm.select(animal_telem[[i]], 
                               CTMM = ctmm_guess_i, 
                               verbose = T)
  ctmm_select_summ <- summary(ctmm_select_i)
  saveRDS(ctmm_select_i, paste("ExampleScripts/Output/ctmms/ctmm_select/", names(animal_telem)[i], ".rds", sep = ""))
  saveRDS(ctmm_select_summ, paste("ExampleScripts/Output/ctmms/ctmm_select_summ/", names(animal_telem)[i], ".rds", sep = ""))
  end <- proc.time() # record the time at which this line executes
  print(paste("Animal ", i, " of ", length(animal_telem),
              " complete. Runtime = ", end - start, ".", sep = ""))
}


## B. Extract estimate and IQR for tau_position, tau_velocity, area, speed ----

selected_model <- tau_position <- tau_velocity <- omega <- tau_position_est <- area_est <- tau_velocity_est <- speed_est <- rep(NA, length(animal_telem))
for(i in 9:length(animal_telem)){ # note: 8 is missing, just skip. 
  k <- readRDS(paste("ExampleScripts/Output/ctmms/ctmm_select/", names(animal_telem)[i], ".rds", sep = ""))
  k_summ <- readRDS(paste("ExampleScripts/Output/ctmms/ctmm_select_summ/", names(animal_telem)[i], ".rds", sep = ""))
  selected_model[i] <- rownames(k_summ)[which.min(k_summ[, 1])]
  k_summ_new <- summary(k[[which.min(k_summ[, 1])]])
  tau_position_est[i] <- k_summ_new$CI[2, 2]
  area_est[i] <- k_summ_new$CI[1, 2]
  tau_velocity_est[i] <- k_summ_new$CI[3, 2]
  speed_est[i] <- k_summ_new$CI[4, 2]
  tau_position[i] <- k[[which.min(k_summ[, 1])]]$tau[1]
  tau_velocity[i] <- k[[which.min(k_summ[, 1])]]$tau[2]
  omega[i] <- k[[which.min(k_summ[, 1])]]$omega
  rm(k)
  rm(k_summ)
}

## C. Summarize outputs across animals ----
best_model <- table(selected_model)
best_model

speed_iqr_kmday <- quantile(na.omit(speed_est), c(0.025, 0.25, 0.5, 0.75, 0.975))
speed_iqr_kmday

area_iqr_kmsq <- quantile(na.omit(area_est), c(0.025, 0.25, 0.5, 0.75, 0.975))
area_iqr_kmsq

tau_posn_iqr_months <- quantile(na.omit(tau_position_est), c(0.025, 0.25, 0.5, 0.75, 0.975))
tau_posn_iqr_months

tau_velo_iqr_months <- quantile(na.omit(tau_velocity_est), c(0.025, 0.25, 0.5, 0.75, 0.975))
tau_velo_iqr_months

tab_summary <- rbind(area_iqr_kmsq, 
                     speed_iqr_kmday, 
                     tau_posn_iqr_months, 
                     tau_velo_iqr_months)
tab_summary
