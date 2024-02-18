## Introduction to Animal Movement Ecology
## K. Manlove, August 2023
##
## Path segmentation
## 
## Aspects of this code were drawn from code developed by Simona Picardi.
## Aspects of this code were drawn from code developed by Lauren Ricci. 

# I. Read in libraries and prep data ----
## A. libraries ----
library(adehabitatLT)
library(bcpa)
library(RColorBrewer)
library(momentuHMM)

## B. data ----
bookcliffs <- read.csv("ExampleScripts/Data/UDWR_Elk_BookCliffs_Cleaned.csv", header = T)

# convert dates to POSIXct
bookcliffs$datePOSIX <- lubridate::ymd_hms(bookcliffs$t_)
bookcliffs$animal_id <- factor(bookcliffs$animal_id)
bookcliffs$latitude <- bookcliffs$latitude # be sure your latitudes are contained in a vector
# named "latitude" with lower-case l
bookcliffs$longitude <- bookcliffs$longitude 

bookcliffs$ID <- bookcliffs$animal_id
bookcliffs$x <- bookcliffs$longitude
bookcliffs$y <- bookcliffs$latitude


elk_in <- levels(factor(bookcliffs$animal_id))

# II. Fit Bayesian changepoint models to the data ----
# set up moveData object
data_in <- bookcliffs
animalIn <- elk_in[1]

## A. getCPs function extracts changepoints ----
getCPs <- function(animalIn,
                   data_in){
  # cut down to just points for this animal
  animal_pts <- subset(data_in, animal_id == animalIn)
  
  # convert points into ltraj -- BE CAREFUL that names match everywhere!
  ltraj_animal <- as.ltraj(xy = animal_pts[, c("longitude", "latitude")], 
                           date = animal_pts$datePOSIX,
                           id = "animal_id", 
                           typeII = TRUE)
  ltraj_sm <- ltraj_animal[[1]]
  # calculate velocity along ltraj
  ltraj_sm$vel <- ltraj_sm$dist / ltraj_sm$dt
  
  # build points for getting the VT
  animal_getVT <- animal_pts
  animal_getVT$X <- animal_getVT$longitude
  animal_getVT$Y <- animal_getVT$latitude
  animal_getVT$Time <- animal_getVT$datePOSIX
  vt_animal <- bcpa::GetVT(Data = animal_getVT, units = "day")
  vt_animal$theta <- vt_animal$Theta
  
  # sweep over timeseries for each individual to find significant changepoints
  windowSw_animal_persVel_run03 <- bcpa::WindowSweep(data = vt_animal,
                                                     variable = "V*cos(theta)",
                                                     windowsize = 20,
                                                     windowstep = 1, 
                                                     K = 2,
                                                     tau = T,
                                                     range = .6,
                                                     progress = T
  )
  
  
  # summarize changepoints
  cp <- bcpa::ChangePointSummary(windowSw_animal_persVel_run03)
  
  # pay attention to the modelmode.
  # 1 = mu changed (higher mu = faster, more directed movt)
  # 2 = sigma changed (higher sigma = more tortuous movt)
  # 3 = rho changed (higher rho = more autocorrelated movt)
  # 4 = mu & sigma
  # 5 = mu & rho
  # 6 = rho, & sigma
  # 7 = mu, rho, sigma
  
  
  cp$phases$colIn <- ifelse(cp$phases$t0 <= june01Julian - 8, "blue",
                            ifelse(cp$phases$t0 <= aug01Julian - 8, "green",
                                   "red"))
  
  # configure output
  outList <- list(cp = cp, 
                  windowSw_animal_persVel_run03 = windowSw_animal_persVel_run03)
  
  # go back to the ltraj, and this time, prep the data for an HMM
  moveDataMaterial <- moveHMM::prepData(subset(ltraj_animal[[1]], 
                                               select = c("dist", "rel.angle", "x", "y")),
                                        type = "LL")
  moveDataMaterial$ID <- rep(as.character(animalIn), nrow(moveDataMaterial))
  names(moveDataMaterial) <- c("step", "angle", "x", "y", "ID")
  
  # return prepped data. 
  return(list(outList = outList,
              moveDataMaterial = moveDataMaterial, 
              ltraj_animal = ltraj_animal))
}

## B. Look at a single animal
animalIn <- elk_in[2]

# III. Fit a two-state HMM using momentuHMM ----
## A. prepare data ----
### i. Vsn 1: regularly space, low-error data using prepData ----
#### extract data for focal animal
animal_pts <- subset(bookcliffs, animal_id == animalIn)
#### convert data to an ltraj object
ltraj_animal <- as.ltraj(xy = animal_pts[, c("longitude", "latitude")],
                         date = animal_pts$datePOSIX,
                         id = "animal_id", typeII = TRUE)
#### reduce down to just the timeseries of fixes
ltraj_sm <- ltraj_animal[[1]]
#### calculate velocity and append to timeseries
ltraj_sm$vel <- ltraj_sm$dist / ltraj_sm$dt
#### build the moveData object using prepData
prepped_data <- momentuHMM::prepData(subset(ltraj_animal[[1]],
                                            select = c("dist", "rel.angle", "x", "y")),
                                     type = "LL") # LL here is for lat-long
#### plot the prepared data.
plot(prepped_data, compact = T)

### ii. Vsn 2: irregularly timed points, or points with relevant measurement error ----
#### In this case, we'll use momentuHMM::crawlWrap() to crawl the data and perform
#### (potentially multiple) imputations.

crawled_data <- crawlWrap(obsData = subset(ltraj_animal[[1]],
                                           select = c("dist", "rel.angle", "x", "y")),
                          Time.name = "date",
                          timeStep = "2 hours",
                          type = "LL") # LL here is for lat-long

animal_pts$ID <- animal_pts$animal_id
animal_pts$x <- animal_pts$longitude
animal_pts$y <- animal_pts$latitude

# project to UTM coordinates using package rgdal
library(rgdal)
llcoord <- SpatialPoints(animal_pts[, 3:4],
                         proj4string=CRS("+proj=longlat +datum=WGS84"))
utmcoord <- spTransform(llcoord,CRS("+proj=utm +zone=12 ellps=WGS84"))
# add UTM locations to data frame
animal_pts$x <- attr(utmcoord,"coords")[,1]
animal_pts$y <- attr(utmcoord,"coords")[,2]

# crawl the data (generate imputations using the continuous-time movement model
# of Johnson et al. 2008)
crawled_data <- momentuHMM::crawlWrap(obsData = animal_pts,
                                      Time.name = "datePOSIX",
                                      timeStep = "2 hours",
                                      type = "LL") # LL here is for lat-long

# inspect the head of the crawled data
head(crawled_data)

## B. create momentuHMMData object from crwData object ----
momentuhmm_data <- momentuHMM::prepData(data = crawled_data)

# check the ACF plot (autocorrelation plot) on the momentuhmm data. 
acf(momentuhmm_data$step[!is.na(momentuhmm_data$step)], lag.max=100)

# add cosinor covariate based on hour of day
momentuhmm_data$hour <- as.integer(strftime(momentuhmm_data$datePOSIX, format = "%H", tz="GMT"))

## C. fit the two-state model ----
### i. label states and set distributions ----
state_names <- c("encamped","exploratory")
## specify the distributions for the observations (step and angle)
dist = list(step = "gamma", angle = "vm") # note: vm is for von Mises

### ii. set initial parameters ----
# we'll specify a mean and sd for the step length, and then
# calculate the gamma parameters for that disribution below. 
# shape and scale:
source("ExampleScripts/05b_PathSegmentationSource.R")
# Extract gamma parameters for short state (mean step length = 10, sd = 5)
# first, examine the data:
hist(momentuhmm_data$step)
# based on the histogram, I want one state with pretty short steps 
# (maybe length around 100m) and another with longer steps (maybe around 750m)
test_short <- gamma_pars(mean = 100, sd = 50)
test_short
# long state
test_long <- gamma_pars(mean = 500, sd = 100)
test_long

# plot proposed distributions of steplengths within states
hist(momentuhmm_data$step, breaks = 100, freq = FALSE, ylim = c(0, .01),
     las = 1, xlab = "step length", ylab = "relative frequency", main = "",
     xlim = c(0, 2000))
lines(dgamma(x = c(0:2000), shape = test_short[2], 
             rate = test_short[1]), col = "red", lwd = 2)
lines(dgamma(x = c(0:2000), shape = test_long[2], 
             rate = test_long[1]), col = "green", lwd = 2)


Par0_m1 <- list(step = c(100, 500, # gamma shape parameters for the two states
                         50, 100, # gamma rate parameters for the two states
                         .01, .01), # gamma zero-mass parameters for the two states.
                # I think there should be very little zero-inflation in the model, so
                # setting those zero-mass parameters very low. 
                angle = c(0.1,0.9)) # 
### iii. fit v1-v3 of the model ----
m1 <- momentuHMM::fitHMM(data = momentuhmm_data, 
                         nbStates = 2, # fit two-state model
                         dist = dist, # use our specified distributions
                         Par0 = Par0_m1, # use our specified input parameters
                         estAngleMean = list(angle=FALSE), # no need to estimate the angle
                         stateNames = state_names) # name the states according to the names we stated above. 

pr1 <- pseudoRes(m1)
plotPR(m1)
acf(pr1$stepRes[!is.na(pr1$stepRes)], lag.max = 100)

# let transition probabilities vary systematically with hour in day 
formula <- ~ cosinor(hour, period = 24)
# initial parameters (obtained from nested model m1)
Par0_m2 <- getPar0(model=m1, formula=formula)
# fit model
m2 <- fitHMM(data = momentuhmm_data, 
             nbStates = 2, 
             dist = dist, 
             Par0 = Par0_m2$Par,
             beta0=Par0_m2$beta, 
             stateNames = state_names, 
             formula = formula)

plot(m2)
pr2 <- pseudoRes(m2)
acf(pr2$stepRes[!is.na(pr2$stepRes)], lag.max = 100)
plotPR(m2)


# formulas for parameters of state-dependent observation distributions
# (DM stands for data matrix)
DM <- list(step = list(mean = ~ cosinor(hour, period = 24),
                       sd = ~ cosinor(hour, period = 24),
                       zeromass = ~ cosinor(hour, period = 24)),
           angle = list(concentration = ~ cosinor(hour, period = 24)))
# initial parameters (obtained from nested model m2)
Par0_m3 <- getPar0(model = m2, formula = formula, DM = DM)
# fit model
m3 <- fitHMM(data = momentuhmm_data, nbStates = 2, 
             dist = dist, 
             Par0 = Par0_m3$Par,
             beta0 = Par0_m3$beta, DM = DM, stateNames = state_names,
             formula = formula)
plotPR(m3)
plot(m3)

### iv. compare models with AIC ----
AIC(m1, m2, m3)

### v. decode most likely state sequence ("label states") ----
# (this essentially "labels" every step as one state of another)
states <- viterbi(m3)
# table the proportion of time the animal spent in each state
table(states)/nrow(momentuhmm_data)

### vi. summarize model visually ----
# generate many plots of the HMM output
plot(m3, plotCI = TRUE, covs = data.frame(hour = 12))
# compute pseudo-residuals for the steps and the angles
pr3 <- pseudoRes(m3)
# plot the ACF of step pseudo-residuals
acf(pr3$stepRes[!is.na(pr3$stepRes)], lag.max = 100)

par(mfrow = c(1,3))
acf(pr1$stepRes[!is.na(pr1$stepRes)], las = 1,
    lag.max = 100, ylim = c(-.1, 1), main = "Residual ACF (Model 1)")
abline(h = .2, lty = 2, col = "red")
acf(pr2$stepRes[!is.na(pr2$stepRes)], las = 1, main = "Residual ACF (Model 2)", 
    lag.max = 100, ylim = c(-.1, 1))
abline(h = .2, lty = 2, col = "red")
acf(pr3$stepRes[!is.na(pr3$stepRes)], las = 1, main = "Residual ACF (Model 3)",
    lag.max = 100, ylim = c(-.1, 1))
abline(h = .2, lty = 2, col = "red")



# IV. Fit the same data using a twelve-hour model ----
# crawl the data (generate imputations using the continuous-time movement model
# of Johnson et al. 2008)
crawled_data_12hr <- momentuHMM::crawlWrap(obsData = animal_pts,
                                           Time.name = "datePOSIX",
                                           timeStep = "12 hours",
                                           type = "LL") # LL here is for lat-long

# inspect the head of the crawled data
head(crawled_data_12hr)

# create momentuHMMData object from crwData object
momentuhmm_data_12hr <- momentuHMM::prepData(data = crawled_data_12hr)

# add cosinor covariate based on hour of day
momentuhmm_data_12hr$month <- as.integer(strftime(momentuhmm_data_12hr$datePOSIX, format = "%m", tz="GMT"))

# check the ACF plot (autocorrelation plot) on the momentuhmm data. 
acf(momentuhmm_data_12hr$step[!is.na(momentuhmm_data_12hr$step)], lag.max=100)

# fit the two-state model
## label states
state_names <- c("encamped","exploratory")
## specify the distributions for the observations (step and angle)
dist = list(step = "gamma", angle = "vm") # note: vm is for von Mises
# initial parameters
# we'll specify a mean and sd for the step length, and then
# calculate the gamma parameters for that disribution below. 
# shape and scale:
source("ExampleScripts/05b_PathSegmentationSource.R")
# Extract gamma parameters for short state (mean step length = 10, sd = 5)
# first, examine the data:
hist(momentuhmm_data_12hr$step)
# based on the histogram, I want one state with pretty short steps 
# (maybe length around 100m) and another with longer steps (maybe around 750m)
test_short <- gamma_pars(mean = 500, sd = 300)
test_short
# long state
test_long <- gamma_pars(mean = 1000, sd = 200)
test_long

# plot proposed distributions of steplengths within states
hist(momentuhmm_data_12hr$step, breaks = 100, freq = FALSE, ylim = c(0, .0025),
     las = 1, xlab = "step length", ylab = "relative frequency", main = "",
     xlim = c(0, 5000))
lines(dgamma(x = c(0:5000), shape = test_short[2], 
             rate = test_short[1]), col = "red", lwd = 2)
lines(dgamma(x = c(0:5000), shape = test_long[2], 
             rate = test_long[1]), col = "green", lwd = 2)

Par0_m1_12hr <- list(step = c(500, 1000, # gamma shape parameters for the two states
                              300, 200 # gamma rate parameters for the two states
), # gamma zero-mass parameters for the two states.
# I think there should be very little zero-inflation in the model, so
# setting those zero-mass parameters very low. 
angle = c(0.1,0.9)) # 
# fit model
m1_12hr <- momentuHMM::fitHMM(data = momentuhmm_data_12hr, 
                              nbStates = 2, # fit two-state model
                              dist = dist, # use our specified distributions
                              Par0 = Par0_m1_12hr, # use our specified input parameters
                              estAngleMean = list(angle=FALSE), # no need to estimate the angle
                              stateNames = state_names) # name the states according to the names we stated above. 

pr1_12hr <- pseudoRes(m1_12hr)
acf(pr1_12hr$stepRes[!is.na(pr1_12hr$stepRes)], lag.max = 100)
plotPR(m1_12hr)
# note how good this ACF looks! We're in a pretty space with respect to
# residual temporal autocorrelation. 

# let transition probabilities vary systematically with hour in day 
formula_12hr <- ~ cosinor(month, period = 12)
# initial parameters (obtained from nested model m1)
Par0_m2_12hr <- getPar0(model=m1_12hr, formula=formula_12hr)
# fit model
m2_12hr <- fitHMM(data = momentuhmm_data_12hr, 
                  nbStates = 2, 
                  dist = dist, 
                  Par0 = Par0_m2_12hr$Par,
                  beta0=Par0_m2_12hr$beta, 
                  stateNames = state_names, 
                  formula = formula)

plotPR(m2_12hr)
plot(m2_12hr)
pr2_12hr <- pseudoRes(m2_12hr)
acf(pr2_12hr$stepRes[!is.na(pr2_12hr$stepRes)], lag.max = 100)


# formulas for parameters of state-dependent observation distributions
# (DM stands for data matrix)
DM_12hr <- list(step = list(mean = ~ cosinor(month, period = 12),
                            sd = ~ cosinor(month, period = 12)),
                angle = list(concentration = ~ cosinor(month, period = 12)))
# initial parameters (obtained from nested model m2)
Par0_m3_12hr <- getPar0(model = m2_12hr, formula = formula_12hr, DM = DM_12hr)
# fit model
m3_12hr <- fitHMM(data = momentuhmm_data_12hr, nbStates = 2, 
                  dist = dist, 
                  Par0 = Par0_m3_12hr$Par,
                  beta0 = Par0_m3_12hr$beta, DM = DM_12hr, stateNames = state_names,
                  formula = formula_12hr)

plot(m3_12hr)

AIC(m1_12hr, m2_12hr, m3_12hr)


# decode most likely state sequence
# (this essentially "labels" every step as one state of another)
states_12hr <- viterbi(m3_12hr)
# table the proportion of time the animal spent in each state
table(states_12hr)/nrow(momentuhmm_data_12hr)
# compute pseudo-residuals for the steps and the angles
pr3_12hr <- pseudoRes(m3_12hr)
# plot the ACF of step pseudo-residuals
acf(pr3_12hr$stepRes[!is.na(pr3_12hr$stepRes)], lag.max = 100)

par(mfrow = c(1,3))
acf(pr1$stepRes[!is.na(pr1$stepRes)], las = 1,
    lag.max = 100, ylim = c(-.1, 1), main = "Residual ACF (Model 1)")
abline(h = .2, lty = 2, col = "red")
acf(pr2$stepRes[!is.na(pr2$stepRes)], las = 1, main = "Residual ACF (Model 2)", 
    lag.max = 100, ylim = c(-.1, 1))
abline(h = .2, lty = 2, col = "red")
acf(pr3$stepRes[!is.na(pr3$stepRes)], las = 1, main = "Residual ACF (Model 3)",
    lag.max = 100, ylim = c(-.1, 1))
abline(h = .2, lty = 2, col = "red")


# V. fit a three-state model ----
## A. set starting parameters for states ----
hist(momentuhmm_data_12hr$step)
# based on the histogram, I want one state with pretty short steps 
# (maybe length around 100m) and another with longer steps (maybe around 750m)
test_short <- gamma_pars(mean = 300, sd = 300)
test_short
# intermediate state
test_inter <- gamma_pars(mean = 600, sd = 200)
test_inter
# long state
test_long <- gamma_pars(mean = 1200, sd = 300)
test_long

# plot proposed distributions of steplengths within states
hist(momentuhmm_data_12hr$step, breaks = 100, freq = FALSE, ylim = c(0, .003),
     las = 1, xlab = "step length", ylab = "relative frequency", main = "",
     xlim = c(0, 5000))
lines(dgamma(x = c(0:5000), shape = test_short[2], 
             rate = test_short[1]), col = "red", lwd = 2)
lines(dgamma(x = c(0:5000), shape = test_inter[2], 
             rate = test_inter[1]), col = "blue", lwd = 2)
lines(dgamma(x = c(0:5000), shape = test_long[2], 
             rate = test_long[1]), col = "green", lwd = 2)


state_names_3state <- c("encamped", "drifting", "exploratory")

## B. Build HMM ----
Par0_m1_12hr_3state <- list(step = c(300, 600, 1200,# gamma shape parameters for the two states
                                     300, 200, 300 # gamma rate parameters for the two states
), # gamma zero-mass parameters for the two states.
# I think there should be very little zero-inflation in the model, so
# setting those zero-mass parameters very low. 
angle = c(0.1, 0.9, 1)) # 
# fit model
m1_12hr_3state <- momentuHMM::fitHMM(data = momentuhmm_data_12hr, 
                                     nbStates = 3, # fit three-state model
                                     dist = dist, # use our specified distributions
                                     Par0 = Par0_m1_12hr_3state, # use our specified input parameters
                                     estAngleMean = list(angle=FALSE), # no need to estimate the angle
                                     stateNames = state_names_3state) # name the states according to the names we stated above. 

## C. examine HMM ----
pr1_12hr_3state <- pseudoRes(m1_12hr_3state)
acf(pr1_12hr_3state$stepRes[!is.na(pr1_12hr_3state$stepRes)], lag.max = 100)
plotPR(m1_12hr_3state)

# ACF looks pretty good...
plot(m1_12hr_3state)
# this actually looks pretty good -- the exploratory state is picking off the really big moves....

## D. Compare to the two-state models ----
AIC(m1_12hr, m2_12hr, m3_12hr, m1_12hr_3state)
# Even this "simple" three-state is the winner in terms of AIC, 
# at least for this animal. 



# VI. Fit and assess separate models for a bunch of animals ----
## A. build function to fit multiple HMMs ----
elkHMM <- function(animalIn, dataIn){
  # prep and crawl animal's data
  animal_pts <- subset(bookcliffs, animal_id == animalIn)
  ## project to UTM coordinates using package rgdal
  llcoord <- SpatialPoints(animal_pts[, 3:4],
                           proj4string=CRS("+proj=longlat +datum=WGS84"))
  utmcoord <- spTransform(llcoord,CRS("+proj=utm +zone=12 ellps=WGS84"))
  ## add UTM locations to data frame
  animal_pts$x <- attr(utmcoord,"coords")[,1]
  animal_pts$y <- attr(utmcoord,"coords")[,2]
  
  # crawl the data
  crawled_data <- momentuHMM::crawlWrap(obsData = animal_pts,
                                        Time.name = "datePOSIX",
                                        timeStep = "12 hours",
                                        type = "LL") # LL here is for lat-long
  
  ## create momentuHMMData object from crwData object ----
  momentuhmm_data <- momentuHMM::prepData(data = crawled_data)
  
  ## Build HMM ----
  Par0_m1_12hr_3state <- list(step = c(300, 600, 1200,# gamma shape parameters for the two states
                                       300, 200, 300 # gamma rate parameters for the two states
  ), # gamma zero-mass parameters for the two states.
  # I think there should be very little zero-inflation in the model, so
  # setting those zero-mass parameters very low. 
  angle = c(0.1, 0.9, 1)) # 
  # fit model
  m1_12hr_3state <- momentuHMM::fitHMM(data = momentuhmm_data_12hr, 
                                       nbStates = 3, # fit three-state model
                                       dist = dist, # use our specified distributions
                                       Par0 = Par0_m1_12hr_3state, # use our specified input parameters
                                       estAngleMean = list(angle=FALSE), # no need to estimate the angle
                                       stateNames = state_names_3state) # name the states according to the names we stated above. 
  
  # aggregate output param estimates:
  m_3state_steps <- m1_12hr_3state$mle$step[1, ]
  m_3state_angles <- m1_12hr_3state$mle$angle[1, ]
  
  # decode steps for output
  steps_decoded <- viterbi(m1_12hr_3state)
  data_with_decoded_steps <- as.data.frame(cbind(momentuhmm_data_12hr, steps_decoded))
  
  # extract MLEs for step length and turning angle in each state
  step_list <- list(m_3state_steps = m_3state_steps)
  angle_list <- list(m_3state_angles = m_3state_angles)
  
  # extract log-likelihoods
  lnLvec <- c(m_3state$mod$minimum)
  
  # plot and store decoded trajectory
  png(paste("ExampleScripts/Output/path_segmentation/all_elk/threestate_", animal_pts$animal_id[1], ".png", sep = ""),
      height = 480, width = 440)
  plot(elk01_fit$m_3state, cumul = F, 
       plotStationary = F, ask = F,
       plotTracks = T)
  dev.off()
  
  outList <- list(animalIn = animalIn,
                  lnLvec = lnLvec,
                  step_list = step_list,
                  angle_list = angle_list,
                  m_3state = m1_12hr_3state,
                  data_with_decoded_steps = data_with_decoded_steps)
  
  return(outList)
}

# test the function
elk01_fit <- elkHMM(animalIn = elk_in[1], dataIn = bookcliffs)
elk01_fit$data_with_decoded_steps

# step up storage objects to store output from each animal (I always handle these as lists)
individ_fit <- steps <- angles <- individ_decoded <- vector("list", length(elk_in))
skips <- c(38, 86, 122) # these three failed on an early try, so now I'm omitting them
to_run <- seq(1:length(elk_in))[-skips] # remove the skips from the vector of animals
# to run.
for(i in to_run[100]:to_run[length(to_run)]){ # loop over run-able animals
  individ_fit[[i]] <- elkHMM(animalIn = elk_in[i], dataIn = bookcliffs) # fit the HMM
  steps[[i]] <- as.vector(individ_fit[[i]]$step_list$m_3state_steps) # pull off steps
  angles[[i]] <- as.vector(individ_fit[[i]]$angle_list$m_3state_angles) # pull off angles
  individ_decoded[[i]] <- individ_fit[[i]]$data_with_decoded_steps # store each individual's decoded trajectory
  print(i)
}


step_summary <- do.call("rbind", steps)
angle_summary <- do.call("rbind", angles)

animal_summary <- as.data.frame(cbind(step_summary, angle_summary))
names(animal_summary) <- c("encamped_step", "intermediate_step", "exploratory_step",
                           "encamped_angle", "intermediate_angle", "exploratory_angle")

decoded_data <- as.data.frame(do.call("rbind", individ_decoded))
dim(decoded_data)

## B. Plots ----
### i. states by animal ----
### momentuHMM isn't using turning angle in its classifier. So, I'm going to
### estimate the von Mises parameters for each state retrospectively. 
encamped_steps <- subset(decoded_data, steps_decoded == 1)
library(CircStats)
n_individ <- length(levels(factor(encamped_steps$ID)))
encamped_angle <- intermed_angle <- explore_angle <- rep(NA, n_individ)
for(i in 1:n_individ){
  individ <- subset(decoded_data, as.character(ID) == levels(factor(decoded_data$ID))[i])
  individ_encamped <- subset(individ, steps_decoded == 1)
  encamped_angle[i] <- vm.ml(na.omit(encamped_steps$angle))
  #encamped_angle_ci <- vm.bootstrap.ci(na.omit(encamped_steps$angle))
}

plot(animal_summary$)

### ii. which states when?



# ### i. two-state model ----
# ## initial parameters for gamma and von Mises distributions
# ## Simona's better way:
# # Source in the function to convert from mean and standard deviation to gamma 
# # shape and scale:
# source("ExampleScripts/05b_PathSegmentationSource.R")
# # Extract gamma parameters for short state (mean step length = 10, sd = 5)
# test1 <- gamma_pars(mean = 10, sd = 5)
# # Intermediate state
# test2 <- gamma_pars(mean = 150, sd = 75) # 150, 200
# # Long state
# test3 <- gamma_pars(mean = 1000, sd = 500)
# 
# hist(mule$step, breaks = 100, freq = FALSE)
# lines(dgamma(x = c(0:7000), shape = test1[2], rate = test1[1]), col = "red")
# lines(dgamma(x = c(0:7000), shape = test2[2], rate = test2[1]), col = "green")
# lines(dgamma(x = c(0:7000), shape = test3[2], rate = test3[1]), col = "blue")
# 
# 
# 
# 
# 
# 
# 
# ## Kezia's old way:
# mu0 <- c(0.1,1) # step mean (two parameters: one for each state)
# sigma0 <- c(0.1,1) # step SD
# zeromass0 <- c(0.1,0.05) # step zero-mass
# stepPar0 <- c(mu0,sigma0,zeromass0)
# angleMean0 <- c(pi,0) # angle mean
# kappa0 <- c(1,1) # angle concentration
# anglePar0 <- c(angleMean0,kappa0)
# 
# m <- fitHMM(data=moveDataMaterial,
#             nbStates=2,
#             stepPar0=stepPar0,
#             anglePar0=anglePar0)
# m
# par(mfrow = c(1, 2))
# plot(m, plotCI = T)
# 
# ### ii. three-state model ----
# ## initial parameters for gamma and von Mises distributions
# mu0 <- c(0.01, .2, 1) # step mean (two parameters: one for each state)
# sigma0 <- c(0.01, .5, 1) # step SD
# zeromass0 <- c(0.1, 0.1, 0.05) # step zero-mass
# stepPar0 <- c(mu0,sigma0,zeromass0)
# angleMean0 <- c(pi, pi/2, 0) # angle mean
# kappa0 <- c(1,1,1) # angle concentration
# anglePar0 <- c(angleMean0,kappa0)
# 
# m_3state <- fitHMM(data=moveDataMaterial,
#                    nbStates=3,
#                    stepPar0=stepPar0,
#                    anglePar0=anglePar0)
# m_3state
# 
# par(mfrow = c(1, 2))
# plot(m_3state, plotCI = T)
# 
# 
# ## 4 states
# mu0 <- c(0.01, .05, .4, .4) # step mean (two parameters: one for each state)
# sigma0 <- c(0.01, .05, .4, .4) # step SD
# zeromass0 <- c(0.01, 0.01, 0, 0) # step zero-mass
# stepPar0 <- c(mu0,sigma0,zeromass0)
# angleMean0 <- c(pi, pi, pi/2, 0) # angle mean
# kappa0 <- c(1, 1, 1, 1) # angle concentration
# anglePar0 <- c(angleMean0,kappa0)
# 
# m_4state <- fitHMM(data=moveDataMaterial,
#                    nbStates=4,
#                    stepPar0=stepPar0,
#                    anglePar0=anglePar0)
# m_4state
# 
# par(mfrow = c(1, 2))
# plot(m_4state, plotCI = T)
# 
# 
# 
# ### i. prep data for moveHMM ----
# animal_pts_12 <- subset(bookcliffs_12, animal_id == animalIn)
# #& Longitude <= -112.7)
# ltraj_animal_12 <- as.ltraj(xy = animal_pts_12[, c("longitude", "latitude")],
#                             date = animal_pts_12$datePOSIX,
#                             id = "animal_id", typeII = TRUE)
# ltraj_sm_12 <- ltraj_animal_12[[1]]
# ltraj_sm_12$vel <- ltraj_sm_12$dist / ltraj_sm_12$dt
# moveDataMaterial_12 <- moveHMM::prepData(subset(ltraj_animal_12[[1]],
#                                                 select = c("dist", "rel.angle", "x", "y")),
#                                          type = "LL")
# 
# plot(moveDataMaterial_12, compact = T)
# 
# ### i. two-state model ----
# ## initial parameters for gamma and von Mises distributions
# mu0_12 <- c(0.2, 2) # step mean (two parameters: one for each state)
# sigma0_12 <- c(0.1, 1) # step SD
# zeromass0_12 <- c(0.1, 0.05) # step zero-mass
# stepPar0_12 <- c(mu0_12, sigma0_12)
# angleMean0_12 <- c(pi, 0) # angle mean
# kappa0_12 <- c(1, 1) # angle concentration
# anglePar0_12 <- c(angleMean0_12, kappa0_12)
# 
# m_12 <- fitHMM(data = moveDataMaterial_12,
#                nbStates = 2,
#                stepPar0 = stepPar0_12,
#                anglePar0 = anglePar0_12)
# m_12
# par(mfrow = c(1, 2))
# plot(m, plotCI = T)
# 
# 
# ### ii. three-state model ----
# mu0_12_three <- c(0.2, 1, 2) # step mean (two parameters: one for each state)
# sigma0_12_three <- c(1, 0, 0) # step SD
# zeromass0_12_three <- c(0.1, 0.05, 0.05) # step zero-mass
# stepPar0_12_three <- c(mu0_12_three, sigma0_12_three)
# angleMean0_12_three <- c(pi, 0, 0) # angle mean
# kappa0_12_three <- c(1, 1, 1) # angle concentration
# anglePar0_12_three <- c(angleMean0_12_three, kappa0_12_three)
# 
# m_12_three <- fitHMM(data = moveDataMaterial_12,
#                      nbStates = 3,
#                      stepPar0 = stepPar0_12_three,
#                      anglePar0 = anglePar0_12_three)
# m_12_three
# par(mfrow = c(1, 2))
# plot(m, plotCI = T)
# 
# 
# 
# # ## 4 states -- reinitialized
# # mu0 <- c(0.01, .05, .25, .1) # step mean (two parameters: one for each state)
# # sigma0 <- c(0.01, .05, .17, .4) # step SD
# # zeromass0 <- c(0.01, 0.01, 0, 0) # step zero-mass
# # stepPar0 <- c(mu0,sigma0,zeromass0)
# # angleMean0 <- c(pi, pi, pi, 0) # angle mean
# # kappa0 <- c(1, 1, 1, 1) # angle concentration
# # anglePar0 <- c(angleMean0,kappa0)
# # 
# # m_4state2 <- fitHMM(data=moveDataMaterial,
# #                    nbStates=4,
# #                    stepPar0=stepPar0,
# #                    anglePar0=anglePar0)
# # m_4state2
# # 
# # par(mfrow = c(2, 2))
# # plot(m_4state2, plotCI = T)
# # 
# ## 5 states -- reinitialized
# mu0 <- c(0.01, .05, .25, .25, 1) # step mean (two parameters: one for each state)
# sigma0 <- c(0.01, .05, .17, .17, .4) # step SD
# zeromass0 <- c(0.01, 0.01, 0, 0, 0) # step zero-mass
# stepPar0 <- c(mu0,sigma0,zeromass0)
# angleMean0 <- c(pi, 4*pi/5, pi, 0, 4*pi/5) # angle mean
# kappa0 <- c(1, 1, 1, 1, 1) # angle concentration
# anglePar0 <- c(angleMean0,kappa0)
# 
# m_5state <- fitHMM(data = moveDataMaterial,
#                    nbStates = 5,
#                    stepPar0 = stepPar0,
#                    anglePar0 = anglePar0)
# m_5state
# # bind on state values decoded through Viterbi algorithm
# m_5state_states <- viterbi(m_5state)
# ltraj_animal_with_states <- as.data.frame(cbind(ltraj_animal[[1]], states))
# class(moveDataMaterial)
# moveDataWithStates <- cbind(moveDataMaterial, states)


# cut down to just points in states 1 and 2, and map those:
stationary_pts <- moveDataWithStates[which(moveDataWithStates$states %in% c(1, 2, 3)), ]
non_stationary_pts <- moveDataWithStates[!(moveDataWithStates$states %in% c(1, 2, 3)), ]
plot(stationary_pts$y ~ stationary_pts$x)
points(non_stationary_pts$y ~ non_stationary_pts$x, col = "grey80")
moveDataWithStates$state_recode <- ifelse(moveDataWithStates$states %in% c(1, 2, 3),
                                          "stationary", "transient")
# Extract centroid for each block of consecutive stationary points.  
# that is the point that gets used for assessing fidelity
# (idea from Harris,Descamps et al., Jo Anim Ecol 2019)
moveDataWithStates$state_recode_shifted <- c(moveDataWithStates$state_recode[-1], NA)
switches <- which (moveDataWithStates$state_recode != moveDataWithStates$state_recode_shifted)
table(table(switches))
new_stationary <- which(moveDataWithStates$state_recode != moveDataWithStates$state_recode_shifted & 
                          moveDataWithStates$state_recode_shifted == "stationary")
end_stationary <- which(moveDataWithStates$state_recode != moveDataWithStates$state_recode_shifted & 
                          moveDataWithStates$state_recode_shifted == "transient")
stationary_bouts <- vector("list", length(new_stationary))
bouts_center <- data.frame("x" = numeric(length(new_stationary)),
                           "y" = numeric(length(new_stationary)))

# 1) extract stationary blocks
for(i in 1:(length(new_stationary) - 1)){
  stationary_bouts[[i]] <- moveDataWithStates[new_stationary[i]:end_stationary[i], ]
  bouts_center$x[i] <- mean(stationary_bouts[[i]]$x)
  bouts_center$y[i] <- mean(stationary_bouts[[i]]$y)
}

col_in <- rainbow(start = .5, end = 1, n = length(new_stationary))
plot(stationary_pts$y ~ stationary_pts$x)
points(non_stationary_pts$y ~ non_stationary_pts$x, col = "grey80")
points(bouts_center$y ~ bouts_center$x, col = col_in, pch = 16)

#
plotStates(m_5state)
plotStationary(m_5state, plotCI=TRUE)

#param estimates:
m_5state$mle


#lnL:
m_5state$mod$minimum


# #-----------
# #-- functionalized vsn:
# #-----------
# multipleHMMs <- function(animalIn, dataIn){
#   
#   # prep data for moveHMM:
#   animal_pts <- subset(dataIn, Animal == animalIn & Longitude <= -112.7 & Longitude >= -113.0)
#   # pull off the first and last 10 fixes
#   animal_pts <- animal_pts[-c(1:10, (dim(animal_pts)[1] - 10):dim(animal_pts)[1]), ]
#   ltraj_animal <- as.ltraj(xy = animal_pts[, c("Longitude", "Latitude")], 
#                            date = animal_pts$datePOSIX,
#                            id = "Animal", typeII = TRUE)
#   ltraj_sm <- ltraj_animal[[1]]
#   ltraj_sm$vel <- ltraj_sm$dist / ltraj_sm$dt
#   moveDataMaterial <- prepData(subset(ltraj_animal[[1]], 
#                                       select = c("dist", "rel.angle", "x", "y")),
#                                type = "LL")
#   
#   
#   ## initial parameters for gamma and von Mises distributions
#   mu0 <- c(0.1,1) # step mean (two parameters: one for each state)
#   sigma0 <- c(0.1,1) # step SD
#   zeromass0 <- c(0.1,0.05) # step zero-mass
#   stepPar0 <- c(mu0,sigma0,zeromass0)
#   angleMean0 <- c(pi,0) # angle mean
#   kappa0 <- c(1,1) # angle concentration
#   anglePar0 <- c(angleMean0,kappa0)
#   
#   m_2state <- fitHMM(data=moveDataMaterial,
#                      nbStates=2,
#                      stepPar0=stepPar0,
#                      anglePar0=anglePar0)
#   print("m_2state complete")
#   
#   ## 3 states
#   mu0 <- c(0.01, .2, 1) # step mean (two parameters: one for each state)
#   sigma0 <- c(0.01, .5, 1) # step SD
#   zeromass0 <- c(0.1, 0.1, 0.05) # step zero-mass
#   stepPar0 <- c(mu0,sigma0,zeromass0)
#   angleMean0 <- c(pi, pi/2, 0) # angle mean
#   kappa0 <- c(1,1,1) # angle concentration
#   anglePar0 <- c(angleMean0,kappa0)
#   
#   m_3state <- fitHMM(data=moveDataMaterial,
#                      nbStates=3,
#                      stepPar0=stepPar0,
#                      anglePar0=anglePar0)
#   print("m_3state complete")
#   
#   
#   ## 4 states -- reinitialized
#   mu0 <- c(0.01, .05, .25, .1) # step mean (two parameters: one for each state)
#   sigma0 <- c(0.01, .05, .17, .4) # step SD
#   zeromass0 <- c(0.01, 0.01, 0, 0) # step zero-mass
#   stepPar0 <- c(mu0,sigma0,zeromass0)
#   angleMean0 <- c(pi, pi, pi, 0) # angle mean
#   kappa0 <- c(1, 1, 1, 1) # angle concentration
#   anglePar0 <- c(angleMean0,kappa0)
#   
#   m_4state <- fitHMM(data=moveDataMaterial,
#                      nbStates=4,
#                      stepPar0=stepPar0,
#                      anglePar0=anglePar0)
#   print("m_4state complete")
#   
#   
#   ## 5 states -- reinitialized
#   mu0 <- c(0.01, .05, .25, .25, 1) # step mean (two parameters: one for each state)
#   sigma0 <- c(0.01, .05, .17, .17, .4) # step SD
#   zeromass0 <- c(0.01, 0.01, 0, 0, 0) # step zero-mass
#   stepPar0 <- c(mu0,sigma0,zeromass0)
#   angleMean0 <- c(pi, 4*pi/5, pi, 0, 4*pi/5) # angle mean
#   kappa0 <- c(1, 1, 1, 1, 1) # angle concentration
#   anglePar0 <- c(angleMean0,kappa0)
#   
#   m_5state <- fitHMM(data = moveDataMaterial,
#                      nbStates = 5,
#                      stepPar0 = stepPar0,
#                      anglePar0 = anglePar0)
#   print("m_5state complete")
#   
#   
#   # aggregate output param estimates:
#   m_2state_steps <- cbind(t(m_2state$mle$stepPar), rep(2, 2))
#   m_2state_angles <- cbind(t(m_2state$mle$anglePar), rep(2, 2))
#   m_3state_steps <- cbind(t(m_3state$mle$stepPar), rep(3, 3))
#   m_3state_angles <- cbind(t(m_3state$mle$anglePar), rep(3, 3))
#   m_4state_steps <- cbind(t(m_4state$mle$stepPar), rep(4, 4))
#   m_4state_angles <- cbind(t(m_4state$mle$anglePar), rep(4, 4))
#   m_5state_steps <- cbind(t(m_5state$mle$stepPar), rep(5, 5))
#   m_5state_angles <- cbind(t(m_5state$mle$anglePar), rep(5, 5))
#   
#   step_list <- list(m_2state_steps = m_2state_steps,
#                     m_3state_steps = m_3state_steps,
#                     m_4state_steps = m_4state_steps,
#                     m_5state_steps = m_5state_steps)
#   
#   angle_list <- list(m_2state_angles = m_2state_angles,
#                      m_3state_angles = m_3state_angles,
#                      m_4state_angles = m_4state_angles,
#                      m_5state_angles = m_5state_angles)
#   
#   #lnL:
#   m_5state$mod$minimum
#   
#   lnLvec <- c(m_2state$mod$minimum,
#               m_3state$mod$minimum,
#               m_4state$mod$minimum,
#               m_5state$mod_minimum)
#   
#   
#   
#   outList <- list(animalIn = animalIn,
#                   lnLvec = lnLvec,
#                   step_list = step_list,
#                   angle_list = angle_list,
#                   m_2state = m_2state,
#                   m_3state = m_3state,
#                   m_4state = m_4state,
#                   m_5state = m_5state)
#   
#   return(outList)
# }
# 
# ramsIn
# ram_39859_2_hmm <- ramHMMs(animalIn = ramsIn[2])
# ram_39859_2_hmm$lnLvec
# plot(ram_39859_2_hmm$m_5state)
# plot(ram_39859_2_hmm$m_4state)
# ram_39860_2_hmm <- ramHMMs(animalIn = ramsIn[4])
# ram_39860_2_hmm$lnLvec
# plot(ram_39860_2_hmm$m_5state)
# plot(ram_39860_2_hmm$m_4state)
# ram_39861_1_hmm <- ramHMMs(animalIn = ramsIn[5])
# ram_39861_1_hmm$lnLvec
# plot(ram_39861_1_hmm$m_5state)
# plot(ram_39861_1_hmm$m_4state)
# 
# # too little data for 39862_1?
# ram_39862_1_hmm <- ramHMMs(animalIn = ramsIn[6])
# 
# ram_39863_1_hmm <- ramHMMs(animalIn = ramsIn[7])
# ram_39863_1_hmm$lnLvec
# plot(ram_39863_1_hmm$m_5state)
# plot(ram_39863_1_hmm$m_4state)
# 
# ram_39864_1_hmm <- ramHMMs(animalIn = ramsIn[8])
# ram_39864_1_hmm$lnLvec
# plot(ram_39864_1_hmm$m_5state)
# plot(ram_39864_1_hmm$m_4state)
# plot(ram_39864_1_hmm$m_2state)
# 
# 
# ram_39865_2_hmm <- ramHMMs(animalIn = ramsIn[10])
# ram_39866_2_hmm <- ramHMMs(animalIn = ramsIn[12])
# ram_39867_1_hmm <- ramHMMs(animalIn = ramsIn[13])
# ram_39868_1_hmm <- ramHMMs(animalIn = ramsIn[14])
# ram_39868_2_hmm <- ramHMMs(animalIn = ramsIn[15])
# ram_39869_1_hmm <- ramHMMs(animalIn = ramsIn[16])
# ram_39870_1_hmm <- ramHMMs(animalIn = ramsIn[17])
# ram_39871_1_hmm <- ramHMMs(animalIn = ramsIn[18])
# ram_39872_1_hmm <- ramHMMs(animalIn = ramsIn[19])
# ram_39872_2_hmm <- ramHMMs(animalIn = ramsIn[20])
# 
# 
# hmmList <- list(ram_39859_2_hmm,
#                 ram_39860_2_hmm,
#                 ram_39861_1_hmm,
#                 #ram_39862_1_hmm,
#                 ram_39863_1_hmm,
#                 ram_39864_1_hmm,
#                 ram_39865_2_hmm,
#                 ram_39866_2_hmm,
#                 ram_39867_1_hmm,
#                 ram_39868_1_hmm,
#                 ram_39868_2_hmm,
#                 ram_39869_1_hmm,
#                 ram_39871_1_hmm,
#                 ram_39872_1_hmm,
#                 ram_39872_2_hmm)
# # extract all coefficients from the 5-state model, and plot. 

meanStepLength <- meanTurningAngle <- stateNumber <- animalID <- vector("list", length(hmmList))
for(i in 1:length(meanStepLength)){
  meanStepLength[[i]] <- hmmList[[i]]$m_5state$mle$stepPar[1, ]
  meanTurningAngle[[i]] <- hmmList[[i]]$m_5state$mle$anglePar[1, ]
  stateNumber[[i]] <- seq(1, 5)
  animalID[[i]] <- rep(as.character(hmmList[[i]]$animalIn), 5)
}

fullMeanStep <- do.call("c", meanStepLength)
fullMeanTurning <- do.call("c", meanTurningAngle)
fullStateNumber <- do.call("c", stateNumber)
fullAnimalID <- do.call("c", animalID)

fullStates <- as.data.frame(cbind(fullMeanStep, 
                                  fullMeanTurning,
                                  fullStateNumber,
                                  fullAnimalID))
fullStates$fullMeanStep <- as.numeric(as.character(fullStates$fullMeanStep))
fullStates$fullMeanTurning <- as.numeric(as.character(fullStates$fullMeanTurning))

colIn <- rainbow(n = 5, start = .6, end = .2)
colIn <- brewer.pal(n = 5, name = "Set1")
colIn <- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99')
par(cex.lab = 1, cex.axis = .8)
plot(fullStates$fullMeanStep ~ abs(fullStates$fullMeanTurning),
     log = "y",  pch = 16,
     col = colIn[fullStates$fullStateNumber],
     las = 1,
     xlab = "|Mean turning angle|",
     ylab = "Mean step length", cex = 1.5)

legText <- c("Static local", 
             "Low movement local", 
             "Moderate movement local", 
             "Moderate movement directed", 
             "Long movement directed")
legend("bottomleft", legText, col = colIn, pch = rep(16, 5),
       pt.cex = rep(1.5, 5), bty = "n")



#-----
#-- V2 with state probs; without simpler models 
#-----
ramHMMsV2 <- function(animalIn){
  
  # prep data for moveHMM:
  animal_pts <- subset(rams, 
                       Animal == animalIn & Longitude <= -112.7 & Longitude >= -113.0)
  # pull off the first and last 10 fixes
  animal_pts <- animal_pts[-c(1:10, (dim(animal_pts)[1] - 10):dim(animal_pts)[1]), ]
  ltraj_animal <- as.ltraj(xy = animal_pts[, c("Longitude", "Latitude")], 
                           date = animal_pts$datePOSIX,
                           id = "Animal", typeII = TRUE)
  ltraj_sm <- ltraj_animal[[1]]
  ltraj_sm$vel <- ltraj_sm$dist / ltraj_sm$dt
  moveDataMaterial <- prepData(subset(ltraj_animal[[1]], 
                                      select = c("dist", "rel.angle", "x", "y")),
                               type = "LL")
  
  ## 5 states -- reinitialized
  mu0 <- c(0.01, .05, .25, .25, 1) # step mean (two parameters: one for each state)
  sigma0 <- c(0.01, .05, .17, .17, .4) # step SD
  zeromass0 <- c(0.01, 0.01, 0, 0, 0) # step zero-mass
  stepPar0 <- c(mu0,sigma0,zeromass0)
  angleMean0 <- c(pi, 4*pi/5, pi, 0, 4*pi/5) # angle mean
  kappa0 <- c(1, 1, 1, 1, 1) # angle concentration
  anglePar0 <- c(angleMean0,kappa0)
  
  m_5state <- fitHMM(data = moveDataMaterial,
                     nbStates = 5,
                     stepPar0 = stepPar0,
                     anglePar0 = anglePar0)
  print("m_5state complete")
  
  #lnL:
  m_5state$mod$minimum
  
  lnLvec <- c(
    m_5state$mod_minimum)
  
  stateProbs_m5state <- stateProbs(m_5state)
  updatedMoveData <- as.data.frame(cbind(animal_pts, 
                                         moveDataMaterial,
                                         stateProbs_m5state)
  )
  stateProbsOut <- apply(stateProbs(m_5state),
                         1, 
                         which.max)
  
  mixingProps <- table(stateProbsOut)/length(stateProbsOut)
  
  outList <- list(animalIn = animalIn,
                  lnLvec = lnLvec,
                  m_5state = m_5state,
                  mixingProps = mixingProps,
                  stateProbsOut = stateProbsOut,
                  updatedMoveData = updatedMoveData)
  
  return(outList)
}


ramsIn
ram_39859_2_hmmV2 <- ramHMMsV2(animalIn = ramsIn[2])
ram_39860_2_hmmV2 <- ramHMMsV2(animalIn = ramsIn[4])
ram_39861_1_hmmV2 <- ramHMMsV2(animalIn = ramsIn[5])

# too little data for 39862_1?
#ram_39862_1_hmmV2 <- ramHMMsV2(animalIn = ramsIn[6])
ram_39863_1_hmmV2 <- ramHMMsV2(animalIn = ramsIn[7])
ram_39864_1_hmmV2 <- ramHMMsV2(animalIn = ramsIn[8])
ram_39865_2_hmmV2 <- ramHMMsV2(animalIn = ramsIn[10])
ram_39866_2_hmmV2 <- ramHMMsV2(animalIn = ramsIn[12])
ram_39867_1_hmmV2 <- ramHMMsV2(animalIn = ramsIn[13])
ram_39868_1_hmmV2 <- ramHMMsV2(animalIn = ramsIn[14])
ram_39868_2_hmmV2 <- ramHMMsV2(animalIn = ramsIn[15])
ram_39869_1_hmmV2 <- ramHMMsV2(animalIn = ramsIn[16])
ram_39870_1_hmmV2 <- ramHMMsV2(animalIn = ramsIn[17])
ram_39871_1_hmmV2 <- ramHMMsV2(animalIn = ramsIn[18])
ram_39872_1_hmmV2 <- ramHMMsV2(animalIn = ramsIn[19])
ram_39872_2_hmmV2 <- ramHMMsV2(animalIn = ramsIn[20])


hmmListV2 <- list(ram_39859_2_hmmV2,
                  ram_39860_2_hmmV2,
                  ram_39861_1_hmmV2,
                  #ram_39862_1_hmm,
                  ram_39863_1_hmmV2,
                  ram_39864_1_hmmV2,
                  ram_39865_2_hmmV2,
                  ram_39866_2_hmmV2,
                  ram_39867_1_hmmV2,
                  ram_39868_1_hmmV2,
                  ram_39868_2_hmmV2,
                  ram_39869_1_hmmV2,
                  ram_39871_1_hmmV2,
                  ram_39872_1_hmmV2,
                  ram_39872_2_hmmV2)

# extract all coefficients from the 5-state model, and plot. 
meanStepLength <- meanTurningAngle <- stateNumber <- animalID <- mixingProps <- vector("list", length(hmmListV2))
individTrajectors <- vector("list", length(hmmListV2))
for(i in 1:length(meanStepLength)){
  meanStepLength[[i]] <- hmmListV2[[i]]$m_5state$mle$stepPar[1, ]
  meanTurningAngle[[i]] <- hmmListV2[[i]]$m_5state$mle$anglePar[1, ]
  stateNumber[[i]] <- seq(1, 5)
  animalID[[i]] <- rep(as.character(hmmListV2[[i]]$animalIn), 5)
  mixingProps[[i]] <- hmmListV2[[i]]$mixingProps
  
  individTrajectors[[i]] <- hmmListV2[[i]]$updatedMoveData
}


fullMeanStep <- do.call("c", meanStepLength)
fullMeanTurning <- do.call("c", meanTurningAngle)
fullStateNumber <- do.call("c", stateNumber)
fullAnimalID <- do.call("c", animalID)
fullMixingProps <- do.call("c", mixingProps)

fullStates <- as.data.frame(cbind(fullMeanStep, 
                                  fullMeanTurning,
                                  fullStateNumber,
                                  fullAnimalID,
                                  fullMixingProps))
fullStates$fullMeanStep <- as.numeric(as.character(fullStates$fullMeanStep))
fullStates$fullMeanTurning <- as.numeric(as.character(fullStates$fullMeanTurning))
fullStates$fullMixingProps <- as.numeric(as.character(fullStates$fullMixingProps))

fullTrajs <- do.call("rbind", individTrajectors)

colIn <- rainbow(n = 5, start = .6, end = .2)
colIn <- brewer.pal(n = 5, name = "Set1")
colIn <- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99')
par(mfrow = c(1, 2), cex.lab = 1, cex.axis = .8)
plot(fullStates$fullMeanStep ~ abs(fullStates$fullMeanTurning),
     log = "y",  pch = 16,
     col = colIn[fullStates$fullStateNumber],
     #cex = fullStates$fullMixingProps * 6,
     cex = 1.5,
     las = 1,
     xlab = "|Mean turning angle|",
     ylab = "Mean step length")
plot(fullStates$fullMeanStep ~ abs(fullStates$fullMeanTurning),
     log = "y",  pch = 16,
     col = colIn[fullStates$fullStateNumber],
     cex = fullStates$fullMixingProps * 6,
     las = 1,
     xlab = "|Mean turning angle|",
     ylab = "Mean step length")

legText <- c("Static local", 
             "Low movement local", 
             "Moderate movement local", 
             "Moderate movement directed", 
             "Long movement directed")
legend("bottomleft", legText, col = colIn, pch = rep(16, 5),
       pt.cex = rep(1.5, 5), bty = "n")


#------
#-- how does coposition of states change through time?

# consider polar plot (radian.plot in plotrix())
# break up individTrajectors into weeks, calculate composition of
# states within each week.

fullTrajs$state1 <- fullTrajs[, 32]
fullTrajs$state2 <- fullTrajs[, 33]
fullTrajs$state3 <- fullTrajs[, 34]
fullTrajs$state4 <- fullTrajs[, 35]
fullTrajs$state5 <- fullTrajs[, 36]
fullTrajs$weekInYr <- floor(fullTrajs$dayInYr/7)
state1Avg <- tapply(fullTrajs$state1,
                    fullTrajs$weekInYr,
                    mean)
state2Avg <- tapply(fullTrajs$state2,
                    fullTrajs$weekInYr,
                    mean)
state3Avg <- tapply(fullTrajs$state3,
                    fullTrajs$weekInYr,
                    mean)
state4Avg <- tapply(fullTrajs$state4,
                    fullTrajs$weekInYr,
                    mean)
state5Avg <- tapply(fullTrajs$state5,
                    fullTrajs$weekInYr,
                    mean)
par(mfrow = c(1, 1))
plot(state1Avg ~ c(1:365), ylim = c(0, .5))
points(state2Avg ~ c(1:365), col = "red")
points(state3Avg ~ c(1:365), col = "blue")
points(state4Avg ~ c(1:365), col = "green")
points(state5Avg ~ c(1:365), col = "orange")

par(mfrow = c(1, 1))
plot(state1Avg ~ c(0:52), ylim = c(0, .5), type = "b")
points(state2Avg ~ c(0:52), col = "red", type = "b")
points(state3Avg ~ c(0:52), col = "blue", type = "b")
points(state4Avg ~ c(0:52), col = "green", type = "b")
points(state5Avg ~ c(0:52), col = "orange", type = "b")

# by ram agegroup
fullTrajs$ageCat <- ifelse(fullTrajs$ageNum <= 3.8, "young", ifelse(fullTrajs$ageNum <= 6.9, "prime",
                                                                    "old"))
fullTrajs$weekAge <- paste(fullTrajs$ageCat, "_", fullTrajs$weekInYr, sep = "")
state1Avg <- tapply(fullTrajs$state1,
                    fullTrajs$weekAge,
                    mean)
state2Avg <- tapply(fullTrajs$state2,
                    fullTrajs$weekAge,
                    mean)
state3Avg <- tapply(fullTrajs$state3,
                    fullTrajs$weekAge,
                    mean)
state4Avg <- tapply(fullTrajs$state4,
                    fullTrajs$weekInYr,
                    mean)
state5Avg <- tapply(fullTrajs$state5,
                    fullTrajs$weekAge,
                    mean)

alphaIn <- .5
colIn <- c(rgb(0, 0, 0, alpha = alphaIn),
           rgb(1, 0, 0, alpha = alphaIn),
           rgb(0, 1, 0, alpha = alphaIn),
           rgb(0, 0, 1, alpha = alphaIn),
           rgb(0, 1, 1, alpha = 1))
colIn <- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99')

par(mfrow = c(1, 3))
plot(state1Avg[105:156] ~ c(1:52), ylim = c(0, .5), type = "b",
     col = colIn[1], las = 1,
     xlab = "Week of year", ylab = "Pr[Animal in State]",
     main = "Under 4 years old")
points(state2Avg[105:156] ~ c(1:52), col = colIn[2], type = "b")
points(state3Avg[105:156] ~ c(1:52), col = colIn[3], type = "b")
points(state4Avg[105:156] ~ c(1:52), col = colIn[4], type = "b")
points(state5Avg[105:156] ~ c(1:52), col = colIn[5], type = "b",
       lwd = 3, pch = 16, cex = 1.5)

plot(state1Avg[53:104] ~ c(1:52), ylim = c(0, .5), type = "b",
     col = colIn[1], las = 1,
     xlab = "Week of year", ylab = "Pr[Animal in State]",
     main = "4 to 6.9 years")
points(state2Avg[53:104] ~ c(1:52), col = colIn[2], type = "b")
points(state3Avg[53:104] ~ c(1:52), col = colIn[3], type = "b")
points(state4Avg[53:104] ~ c(1:52), col = colIn[4], type = "b")
points(state5Avg[53:104] ~ c(1:52), col = colIn[5], type = "b",
       lwd = 3, pch = 16, cex = 1.5)

plot(state1Avg[1:52] ~ c(1:52), ylim = c(0, .5), type = "b",
     col = colIn[1], las = 1,
     xlab = "Week of year", ylab = "Pr[Animal in State]",
     main = "Over 6.9 years")
points(state2Avg[1:52] ~ c(1:52), col = colIn[2], type = "b")
points(state3Avg[1:52] ~ c(1:52), col = colIn[3], type = "b")
points(state4Avg[1:52] ~ c(1:52), col = colIn[4], type = "b")
points(state5Avg[1:52] ~ c(1:52), col = colIn[5], type = "b",
       lwd = 3, pch = 16, cex = 1.5)




#----
#-- now, for ewes
#----

#-----
#-- V2 with state probs; without simpler models 
#-----
eweHMMs_4state <- function(animalIn, dataIn){
  
  # prep data for moveHMM:
  animal_pts <- subset(dataIn, as.character(Animal) == as.character(animalIn) & Longitude <= -112.7 & Longitude >= -113.0)
  # pull off the first and last 10 fixes
  animal_pts <- animal_pts[-c(1:10, (dim(animal_pts)[1] - 10):dim(animal_pts)[1]), ]
  ltraj_animal <- as.ltraj(xy = animal_pts[, c("Longitude", "Latitude")], 
                           date = animal_pts$datePOSIX,
                           id = "Animal", typeII = TRUE)
  ltraj_sm <- ltraj_animal[[1]]
  ltraj_sm$vel <- ltraj_sm$dist / ltraj_sm$dt
  moveDataMaterial <- prepData(subset(ltraj_animal[[1]], 
                                      select = c("dist", "rel.angle", "x", "y")),
                               type = "LL")
  
  ## 4 states -- reinitialized
  mu0 <- c(0.01, .05, .25, .1) # step mean (two parameters: one for each state)
  sigma0 <- c(0.01, .05, .17, .4) # step SD
  zeromass0 <- c(0.01, 0.01, 0, 0) # step zero-mass
  stepPar0 <- c(mu0,sigma0,zeromass0)
  angleMean0 <- c(pi, pi, pi, 0) # angle mean
  kappa0 <- c(1, 1, 1, 1) # angle concentration
  anglePar0 <- c(angleMean0,kappa0)
  
  m_4state <- fitHMM(data=moveDataMaterial,
                     nbStates=4,
                     stepPar0=stepPar0,
                     anglePar0=anglePar0)
  print("m_4state complete")
  
  #lnL:
  m_4state$mod$minimum
  
  lnLvec <- c(
    m_4state$mod_minimum)
  
  stateProbs_m4state <- stateProbs(m_4state)
  updatedMoveData <- as.data.frame(cbind(animal_pts, 
                                         moveDataMaterial,
                                         stateProbs_m4state)
  )
  stateProbsOut <- apply(stateProbs(m_4state),
                         1, 
                         which.max)
  
  mixingProps <- table(stateProbsOut)/length(stateProbsOut)
  
  outList <- list(animalIn = animalIn,
                  lnLvec = lnLvec,
                  m_4state = m_4state,
                  mixingProps = mixingProps,
                  stateProbsOut = stateProbsOut,
                  updatedMoveData = updatedMoveData)
  
  return(outList)
}

# # 35480 (currently failing)
# znp_035480 <- read.csv("Data/ZionEweLocations_12Nov2019/ZNP_035480.csv", 
#                        header = T)
# znp_035480$Animal <- as.character(znp_035480$COLLARSERIALNUM)
# znp_035480$datePOSIX <- as.POSIXct(strptime(as.character(znp_035480$DateYearAndJulian), 
#                                             format = "%m/%d/%Y, %I:%M %p"))
# da <- as.character(znp_035480$datePOSIX) # Identifies duplicate fixes
# znp_035480 <- znp_035480[-which(duplicated(paste(da, znp_035480$Animal)) == 1), ] # removes duplicate fixes


# 39879 
znp_039879 <- read.csv("Data/ZionEweLocations_12Nov2019/ZNP_039879_clean.csv", 
                       header = T)
znp_039879$Animal <- paste(as.character(znp_039879$COLLARSERIALNUM), "_", as.character(znp_039879$Deployment), sep = "")
# 039879_1 working
znp_039879_1 <- subset(znp_039879, Animal == "39879_1")
znp_039879_1$datePOSIX <- as.POSIXct(strptime(as.character(znp_039879_1$DateYearAndJulian), 
                                              format = "%m/%d/%Y, %I:%M %p"))
#da <- as.character(znp_039879_1$datePOSIX) # Identifies duplicate fixes
#znp_039879_1b <- znp_039879_1[-which(duplicated(paste(da, znp_039879_1$Animal)) == TRUE), ] # removes duplicate fixes

znp_039879_2 <- subset(znp_039879, Animal == "39879_2")
znp_039879_2$datePOSIX <- as.POSIXct(strptime(as.character(znp_039879_2$DateYearAndJulian), 
                                              format = "%m/%d/%Y, %I:%M %p"))
#da <- as.character(znp_039879_2$datePOSIX) # Identifies duplicate fixes
#znp_039879_2 <- znp_039879_2[-which(duplicated(paste(da, znp_039879_2$Animal)) == 1), ] # removes duplicate fixes


# 39880 (working)
znp_039880 <- read.csv("Data/ZionEweLocations_12Nov2019/ZNP_039880_clean.csv", 
                       header = T)
znp_039880$Animal <- as.character(znp_039880$COLLARSERIALNUM)
znp_039880$datePOSIX <- as.POSIXct(strptime(as.character(znp_039880$DateYearAndJulian), 
                                            format = "%m/%d/%Y, %I:%M %p"))
da <- as.character(znp_039880$datePOSIX) # Identifies duplicate fixes
znp_039880 <- znp_039880[-which(duplicated(paste(da, znp_039880$Animal)) == 1), ] # removes duplicate fixes


# 39881 (untested)
znp_039881 <- read.csv("Data/ZionEweLocations_12Nov2019/ZNP_039881_clean.csv", 
                       header = T)
table(znp_039881$Deployment)
znp_039881$Animal <- as.character(znp_039881$COLLARSERIALNUM)
znp_039881$datePOSIX <- as.POSIXct(strptime(as.character(znp_039881$DateYearAndJulian), 
                                            format = "%m/%d/%Y, %I:%M %p"))
da <- as.character(znp_039881$datePOSIX) # Identifies duplicate fixes
znp_039881 <- znp_039881[-which(duplicated(paste(da, znp_039881$Animal)) == 1), ] # removes duplicate fixes


# 39882 (untested)
znp_039882 <- read.csv("Data/ZionEweLocations_12Nov2019/ZNP_039882_clean.csv", 
                       header = T)
table(znp_039882$Deployment)
znp_039882$Animal <- as.character(znp_039882$COLLARSERIALNUM)
znp_039882$datePOSIX <- as.POSIXct(strptime(as.character(znp_039882$DateYearAndJulian), 
                                            format = "%m/%d/%Y, %I:%M %p"))
#da <- as.character(znp_039882$datePOSIX) # Identifies duplicate fixes
#znp_039882 <- znp_039882[-which(duplicated(paste(da, znp_039882$Animal)) == 1), ] # removes duplicate fixes


# 39883 (untested)
znp_039883 <- read.csv("Data/ZionEweLocations_12Nov2019/ZNP_039883_clean.csv", 
                       header = T)
table(znp_039883$Deployment)
znp_039883$Animal <- as.character(znp_039883$COLLARSERIALNUM)
znp_039883$datePOSIX <- as.POSIXct(strptime(as.character(znp_039883$DateYearAndJulian), 
                                            format = "%m/%d/%Y, %I:%M %p"))
da <- as.character(znp_039883$datePOSIX) # Identifies duplicate fixes
znp_039883 <- znp_039883[-which(duplicated(paste(da, znp_039883$Animal)) == 1), ] # removes duplicate fixes


# 39884 (untested)
znp_039884 <- read.csv("Data/ZionEweLocations_12Nov2019/ZNP_039884_clean.csv", 
                       header = T)
table(znp_039884$Deployment)
znp_039884$Animal <- as.character(znp_039884$COLLARSERIALNUM)
znp_039884$datePOSIX <- as.POSIXct(strptime(as.character(znp_039884$DateYearAndJulian), 
                                            format = "%m/%d/%Y, %I:%M %p"))
da <- as.character(znp_039884$datePOSIX) # Identifies duplicate fixes
znp_039884 <- znp_039884[-which(duplicated(paste(da, znp_039884$Animal)) == 1), ] # removes duplicate fixes


# 39885 (untested)
znp_039885 <- read.csv("Data/ZionEweLocations_12Nov2019/ZNP_039885_clean.csv", 
                       header = T)
table(znp_039885$Deployment)
znp_039885$Animal <- as.character(znp_039885$COLLARSERIALNUM)
znp_039885$datePOSIX <- as.POSIXct(strptime(as.character(znp_039885$DateYearAndJulian), 
                                            format = "%m/%d/%Y, %I:%M %p"))
da <- as.character(znp_039885$datePOSIX) # Identifies duplicate fixes
znp_039885 <- znp_039885[-which(duplicated(paste(da, znp_039885$Animal)) == 1), ] # removes duplicate fixes


# 39886 (untested)
znp_039886 <- read.csv("Data/ZionEweLocations_12Nov2019/ZNP_039886_clean.csv", 
                       header = T)
table(znp_039886$Deployment)
znp_039886$Animal <- paste(as.character(znp_039886$COLLARSERIALNUM), "_", 
                           as.character(znp_039886$Deployment), sep = "")
# 0039886_1 working
znp_039886_1 <- subset(znp_039886, Animal == "39886_1")
znp_039886_1$datePOSIX <- as.POSIXct(strptime(as.character(znp_039886_1$DateYearAndJulian), 
                                              format = "%m/%d/%Y, %I:%M %p"))
#da <- as.character(znp_039879_1$datePOSIX) # Identifies duplicate fixes
#znp_039879_1b <- znp_039879_1[-which(duplicated(paste(da, znp_039879_1$Animal)) == TRUE), ] # removes duplicate fixes

znp_039886_2 <- subset(znp_039886, Animal == "39886_2")
znp_039886_2$datePOSIX <- as.POSIXct(strptime(as.character(znp_039886_2$DateYearAndJulian), 
                                              format = "%m/%d/%Y, %I:%M %p"))
#da <- as.character(znp_039879_2$datePOSIX) # Identifies duplicate fixes
#znp_039879_2 <- znp_039879_2[-which(duplicated(paste(da, znp_039879_2$Animal)) == 1), ] # removes duplicate fixes

# znp_039886$Animal <- as.character(znp_039886$COLLARSERIALNUM)
# znp_039886$datePOSIX <- as.POSIXct(strptime(as.character(znp_039886$DateYearAndJulian), 
#                                             format = "%m/%d/%Y, %I:%M %p"))
# da <- as.character(znp_039886$datePOSIX) # Identifies duplicate fixes
# znp_039886 <- znp_039886[-which(duplicated(paste(da, znp_039886$Animal)) == 1), ] # removes duplicate fixes


# 39887 (untested)
znp_039887 <- read.csv("Data/ZionEweLocations_12Nov2019/ZNP_039887_clean.csv", 
                       header = T)
table(znp_039887$Deployment)
znp_039887$Animal <- as.character(znp_039887$COLLARSERIALNUM)
znp_039887$datePOSIX <- as.POSIXct(strptime(as.character(znp_039887$DateYearAndJulian), 
                                            format = "%m/%d/%Y, %I:%M %p"))
da <- as.character(znp_039887$datePOSIX) # Identifies duplicate fixes
znp_039887 <- znp_039887[-which(duplicated(paste(da, znp_039887$Animal)) == 1), ] # removes duplicate fixes


# 41992 
znp_041992 <- read.csv("Data/ZionEweLocations_12Nov2019/ZNP_041992_clean.csv", 
                       header = T)
table(znp_041992$Deployment)
znp_041992$Animal <- as.character(znp_041992$COLLARSERIALNUM)
znp_041992$datePOSIX <- as.POSIXct(strptime(as.character(znp_041992$DateYearAndJulian), 
                                            format = "%m/%d/%Y, %I:%M %p"))
da <- as.character(znp_041992$datePOSIX) # Identifies duplicate fixes
znp_041992 <- znp_041992[-which(duplicated(paste(da, znp_041992$Animal)) == 1), ] # removes duplicate fixes


# 41993 
znp_041993 <- read.csv("Data/ZionEweLocations_12Nov2019/ZNP_041993_clean.csv", 
                       header = T)
table(znp_041993$Deployment)
znp_041993$Animal <- as.character(znp_041993$COLLARSERIALNUM)
posixDates <- vector("list", dim(znp_041993)[1])
znp_041993$datePOSIX <- as.POSIXct(strptime(as.character(znp_041993$DateYearAndJulian), 
                                            format = "%m/%d/%Y, %I:%M %p"))
da <- as.character(znp_041993$datePOSIX) # Identifies duplicate fixes
znp_041993 <- znp_041993[-which(duplicated(paste(da, znp_041993$Animal)) == 1), ] # removes duplicate fixes


# 41994 currently failing. 
znp_041994 <- read.csv("Data/ZionEweLocations_12Nov2019/ZNP_041994_clean.csv", 
                       header = T)
table(znp_041994$Deployment)
znp_041994$Animal <- as.character(znp_041994$COLLARSERIALNUM)
posixDates <- vector("list", dim(znp_041994)[1])
znp_041994$datePOSIX <- as.POSIXct(strptime(as.character(znp_041994$DateYearAndJulian), 
                                            format = "%m/%d/%Y, %I:%M %p"))
da <- as.character(znp_041994$datePOSIX) # Identifies duplicate fixes
znp_041994 <- znp_041994[-which(duplicated(paste(da, znp_041994$Animal)) == 1), ] # removes duplicate fixes


# 41995 (working)
znp_041995 <- read.csv("Data/ZionEweLocations_12Nov2019/ZNP_041995_clean.csv", 
                       header = T)
table(znp_041995$Deployment)
znp_041995$Animal <- as.character(znp_041995$COLLARSERIALNUM)
posixDates <- vector("list", dim(znp_041995)[1])
znp_041995$datePOSIX <- as.POSIXct(strptime(as.character(znp_041995$DateYearAndJulian), 
                                            format = "%m/%d/%Y, %I:%M %p"))
#da <- as.character(znp_041995$datePOSIX) # Identifies duplicate fixes
#znp_041995 <- znp_041995[-which(duplicated(paste(da, znp_041995$Animal)) == 1), ] # removes duplicate fixes



# fit ewe HMMs
# ewe_znp_035480_hmm <- multipleHMMs(animalIn = levels(factor(znp_035480$Animal))[1],
#                                    dataIn = znp_035480)
ewe_znp_039879_1_hmm <- multipleHMMs(animalIn = levels(factor(znp_039879_1$Animal))[1],
                                     dataIn = znp_039879_1)
ewe_znp_039879_2_hmm <- multipleHMMs(animalIn = levels(factor(znp_039879_2$Animal))[1],
                                     dataIn = znp_039879_2)
ewe_znp_039880_hmm <- multipleHMMs(animalIn = levels(factor(znp_039880$Animal))[1],
                                   dataIn = znp_039880)
ewe_znp_041992_hmm <- multipleHMMs(animalIn = levels(factor(znp_041992$Animal))[1],
                                   dataIn = znp_041992)
ewe_znp_041993_hmm <- multipleHMMs(animalIn = levels(factor(znp_041993$Animal))[1],
                                   dataIn = znp_041993)
ewe_znp_041994_hmm <- multipleHMMs(animalIn = levels(factor(znp_041994$Animal))[1],
                                   dataIn = znp_041994)
ewe_znp_041995_hmm <- multipleHMMs(animalIn = levels(factor(znp_041995$Animal))[1],
                                   dataIn = znp_041995)


# revised fit: 4-state only
#ewe_znp_035480_4st <- eweHMMs_4state(animalIn = levels(factor(znp_035480$Animal))[1],
#                                     dataIn = znp_035480)
ewe_znp_039879_1_4st <- eweHMMs_4state(animalIn = levels(factor(znp_039879_1$Animal))[1],
                                       dataIn = znp_039879_1)
ewe_znp_039879_2_4st <- eweHMMs_4state(animalIn = levels(factor(znp_039879_2$Animal))[1],
                                       dataIn = znp_039879_2)
ewe_znp_039880_4st <- eweHMMs_4state(animalIn = levels(factor(znp_039880$Animal))[1],
                                     dataIn = znp_039880)
ewe_znp_039881_4st <- eweHMMs_4state(animalIn = levels(factor(znp_039881$Animal))[1],
                                     dataIn = znp_039881)
ewe_znp_039882_4st <- eweHMMs_4state(animalIn = levels(factor(znp_039882$Animal))[1],
                                     dataIn = znp_039882)
ewe_znp_039883_4st <- eweHMMs_4state(animalIn = levels(factor(znp_039883$Animal))[1],
                                     dataIn = znp_039883)
ewe_znp_039884_4st <- eweHMMs_4state(animalIn = levels(factor(znp_039884$Animal))[1],
                                     dataIn = znp_039884)
ewe_znp_039885_4st <- eweHMMs_4state(animalIn = levels(factor(znp_039885$Animal))[1],
                                     dataIn = znp_039885)
ewe_znp_039886_1_4st <- eweHMMs_4state(animalIn = levels(factor(znp_039886_1$Animal))[1],
                                       dataIn = znp_039886_1)
ewe_znp_039886_2_4st <- eweHMMs_4state(animalIn = levels(factor(znp_039886_2$Animal))[1],
                                       dataIn = znp_039886_2)
ewe_znp_039887_4st <- eweHMMs_4state(animalIn = levels(factor(znp_039887$Animal))[1],
                                     dataIn = znp_039887)
ewe_znp_041992_4st <- eweHMMs_4state(animalIn = levels(factor(znp_041992$Animal))[1],
                                     dataIn = znp_041992)
ewe_znp_041993_4st <- eweHMMs_4state(animalIn = levels(factor(znp_041993$Animal))[1],
                                     dataIn = znp_041993)
ewe_znp_041994_4st <- eweHMMs_4state(animalIn = levels(factor(znp_041994$Animal))[1],
                                     dataIn = znp_041994)
ewe_znp_041995_4st <- eweHMMs_4state(animalIn = levels(factor(znp_041995$Animal))[1],
                                     dataIn = znp_041995)

# examine output
ewe_znp_041992_hmm$lnLvec
plot(ewe_znp_041992_hmm$m_5state)
plot(ewe_znp_041992_hmm$m_4state)

ewe_znp_041993_hmm$lnLvec
plot(ewe_znp_041993_hmm$m_5state)
plot(ewe_znp_041993_hmm$m_4state)

ewe_znp_041995_hmm$lnLvec
plot(ewe_041995_hmm$m_5state)
plot(ewe_041995_hmm$m_4state)

eweHmmList <- list(
  #ewe_znp_035480_4st,
  ewe_znp_039879_1_4st,
  ewe_znp_039879_2_4st,
  ewe_znp_039880_4st,
  ewe_znp_039881_4st,
  ewe_znp_039882_4st,
  ewe_znp_039883_4st,
  ewe_znp_039884_4st,
  ewe_znp_039885_4st,
  ewe_znp_039886_4st,
  ewe_znp_039887_4st,
  ewe_znp_041992_4st,
  ewe_znp_041993_4st,
  ewe_znp_041994_4st,
  ewe_znp_041995_4st
)

# dput(eweHmmList, "eweHmmList_20191113")
eweHmmList <- dget("eweHmmList_20191113")

#-- relabel states to be common across all animals
stateOrderList <- vector("list", length(eweHmmList))
for(i in 1:length(eweHmmList)){
  stateOrderList[[i]] <- ifelse(abs(eweHmmList[[i]]$m_4state$mle$anglePar[1, ]) <= 1.0, 
                                "ModerateDirected",
                                ifelse(eweHmmList[[i]]$m_4state$mle$stepPar[1, ] <= 0.02, 
                                       "StaticLocal",
                                       ifelse(eweHmmList[[i]]$m_4state$mle$stepPar[1, ] <= 0.2,
                                              "LowMovementStatic",
                                              ifelse(abs(eweHmmList[[i]]$m_4state$mle$anglePar[1, ])  > 1.0,
                                                     "ModerateLocal", NA))))
}

#-- structure of steps and angles in states
# extract all coefficients from the 4-state model, and plot. 
meanStepLength <- meanTurningAngle <- stateNumber <- animalID <- mixingProps <- vector("list", length(eweHmmList))
individTrajectors <- vector("list", length(eweHmmList))
indicesIn <- c(1:length(meanStepLength))
for(i in indicesIn){
  meanStepLength[[i]] <- eweHmmList[[i]]$m_4state$mle$stepPar[1, ]
  meanTurningAngle[[i]] <- eweHmmList[[i]]$m_4state$mle$anglePar[1, ]
  stateNumber[[i]] <- seq(1, 4)
  animalID[[i]] <- rep(as.character(eweHmmList[[i]]$animalIn), 4)
  mixingProps[[i]] <- eweHmmList[[i]]$mixingProps
  individTrajectors[[i]] <- subset(eweHmmList[[i]]$updatedMoveData,
                                   select = names(eweHmmList[[9]]$updatedMoveData))
  
  stateVectors <- cbind(individTrajectors[[i]]$`1`, 
                        individTrajectors[[i]]$`2`,
                        individTrajectors[[i]]$`3`,
                        individTrajectors[[i]]$`4`)
  numberModerateDirected <- length(which(stateOrderList[[i]] == "ModerateDirected"))
  numberModerateLocal <- length(which(stateOrderList[[i]] == "ModerateLocal"))
  numberStaticLocal <- length(which(stateOrderList[[i]] == "StaticLocal"))
  numberLowMovtStatic <- length(which(stateOrderList[[i]] == "LowMovementStatic"))
  if (numberModerateDirected == 1){
    individTrajectors[[i]]$ModerateDirected <- stateVectors[, which(stateOrderList[[i]] == "ModerateDirected")]
  } else if(numberModerateDirected >= 2){
    individTrajectors[[i]]$ModerateDirected <- rowSums(stateVectors[ , 
                                                                     which(stateOrderList[[i]] == "ModerateDirected")])
  } else if (numberModerateDirected == 0){
    individTrajectors[[i]]$ModerateDirected <- rep(0, dim(individTrajectors[[i]])[1])
  }
  
  if (numberModerateLocal == 1){
    individTrajectors[[i]]$ModerateLocal <- stateVectors[, which(stateOrderList[[i]] == "ModerateLocal")]
  } else if (numberModerateLocal >= 2){
    individTrajectors[[i]]$ModerateLocal <- rowSums(stateVectors[ ,which(stateOrderList[[i]] == "ModerateLocal")])
  } else if (numberModerateLocal == 0){
    individTrajectors[[i]]$ModerateLocal <- rep(0, dim(individTrajectors[[i]])[1])
  }
  
  if (numberStaticLocal == 1){
    individTrajectors[[i]]$StaticLocal <- stateVectors[, which(stateOrderList[[i]] == "StaticLocal")]
  } else if(numberStaticLocal >= 2){
    individTrajectors[[i]]$StaticLocal <- rowSums(stateVectors[ ,which(stateOrderList[[i]] == "StaticLocal")])
  } else if (numberStaticLocal == 0){
    individTrajectors[[i]]$StaticLocal <- rep(0, dim(individTrajectors[[i]])[1])
  }
  
  if (numberLowMovtStatic == 1){
    individTrajectors[[i]]$LowMovementStatic <- stateVectors[, which(stateOrderList[[i]] == "LowMovementStatic")]
  } else if(numberLowMovtStatic >= 2){
    individTrajectors[[i]]$LowMovementStatic <- rowSums(stateVectors[ ,which(stateOrderList[[i]] == "LowMovementStatic")])
  } else if (numberLowMovtStatic == 0){
    individTrajectors[[i]]$LowMovementStatic <- rep(0, dim(individTrajectors[[i]])[1])
  }
  
  # individTrajectors[[i]]$moderateDirected <- which(stateOrderList[[i]] == "ModerateDirected")
  
}



fullMeanStep <- do.call("c", meanStepLength)
fullMeanTurning <- do.call("c", meanTurningAngle)
fullStateNumber <- do.call("c", stateNumber)
fullAnimalID <- do.call("c", animalID)
fullMixingProps <- do.call("c", mixingProps)

fullStates <- as.data.frame(cbind(fullMeanStep, 
                                  fullMeanTurning,
                                  fullStateNumber,
                                  fullAnimalID,
                                  fullMixingProps))
fullStates$fullMeanStep <- as.numeric(as.character(fullStates$fullMeanStep))
fullStates$fullMeanTurning <- as.numeric(as.character(fullStates$fullMeanTurning))
fullStates$fullMixingProps <- as.numeric(as.character(fullStates$fullMixingProps))

fullTrajs <- do.call("rbind", individTrajectors)


colIn <- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99')
par(mfrow = c(1, 2), cex.lab = 1, cex.axis = .8)
plot(fullStates$fullMeanStep ~ abs(fullStates$fullMeanTurning),
     log = "y",  pch = 16, ylim = c(0.005, 5),
     col = colIn[fullStates$fullStateNumber],
     #cex = fullStates$fullMixingProps * 6,
     cex = 1.5,
     las = 1,
     xlab = "|Mean turning angle|",
     ylab = "Mean step length")
plot(fullStates$fullMeanStep ~ abs(fullStates$fullMeanTurning),
     log = "y",  pch = 16, ylim = c(0.005, 5),
     col = colIn[fullStates$fullStateNumber],
     cex = fullStates$fullMixingProps * 6,
     las = 1,
     xlab = "|Mean turning angle|",
     ylab = "Mean step length")

legText <- c("Static local", 
             "Low movement local", 
             "Moderate movement local", 
             "Moderate movement directed", 
             "Long movement directed")
legend("bottomleft", legText, col = colIn, pch = rep(16, 5),
       pt.cex = rep(1.5, 5), bty = "n")


#------
#-- how does coposition of states change through time?

# consider polar plot (radian.plot in plotrix())
# break up individTrajectors into weeks, calculate composition of
# states within each week.

fullTrajs$dayInYr <- as.numeric(as.character(strftime(paste(as.character(fullTrajs$date), "", 
                                                            as.character(fullTrajs$time), sep = ""), 
                                                      format = "%j")))
fullTrajs$state1 <- fullTrajs[, 40]
fullTrajs$state2 <- fullTrajs[, 41]
fullTrajs$state3 <- fullTrajs[, 42]
fullTrajs$state4 <- fullTrajs[, 43]
#fullTrajs$state5 <- fullTrajs[, 36]

fullTrajs$weekInYr <- floor(fullTrajs$dayInYr/7)

state2Avg <- tapply(fullTrajs$state2,
                    fullTrajs$weekInYr,
                    mean)
state3Avg <- tapply(fullTrajs$state3,
                    fullTrajs$weekInYr,
                    mean)
state4Avg <- tapply(fullTrajs$state4,
                    fullTrajs$weekInYr,
                    mean)


# plot(state1Avg ~ c(0:52), ylim = c(0, .5), type = "b", col = colIn[1])
# points(state2Avg ~ c(0:52), col = colIn[2], type = "b")
# points(state3Avg ~ c(0:52), col = colIn[3], type = "b")
# points(state4Avg ~ c(0:52), col = colIn[4], type = "b")

# by ewe agegroup
fullTrajs$ageCat <- ifelse(fullTrajs$Age <= 3.8, "young", 
                           ifelse(fullTrajs$Age <= 6.9, "prime",
                                  "old"))
fullTrajs$weekAge <- paste(fullTrajs$ageCat, "_", formatC(fullTrajs$weekInYr, width = 2), sep = "")
state1Avg <- tapply(fullTrajs$state1,
                    fullTrajs$weekAge,
                    mean)
state2Avg <- tapply(fullTrajs$state2,
                    fullTrajs$weekAge,
                    mean)
state3Avg <- tapply(fullTrajs$state3,
                    fullTrajs$weekAge,
                    mean)
state4Avg <- tapply(fullTrajs$state4,
                    fullTrajs$weekAge,
                    mean)

alphaIn <- .5
colIn <- c(rgb(0, 0, 0, alpha = alphaIn),
           rgb(1, 0, 0, alpha = alphaIn),
           rgb(0, 1, 0, alpha = alphaIn),
           rgb(0, 0, 1, alpha = alphaIn),
           rgb(0, 1, 1, alpha = 1))
colIn <- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99')

par(mfrow = c(1, 3))
plot(state1Avg[160:212] ~ c(1:53), ylim = c(0, .6), type = "b",
     col = colIn[1], las = 1,
     xlab = "Week of year", ylab = "Pr[Animal in State]",
     main = "Over 6.9 years")
points(state2Avg[160:212] ~ c(1:53), col = colIn[2], type = "b")
points(state3Avg[160:212] ~ c(1:53), col = colIn[3], type = "b")
points(state4Avg[160:212] ~ c(1:53), col = colIn[4], type = "b")


plot(state1Avg[107:159] ~ c(1:53), ylim = c(0, .6), type = "b",
     col = colIn[1], las = 1,
     xlab = "Week of year", ylab = "Pr[Animal in State]",
     main = "Under 4 years old")
points(state2Avg[107:159] ~ c(1:53), col = colIn[2], type = "b")
points(state3Avg[107:159] ~ c(1:53), col = colIn[3], type = "b")
points(state4Avg[107:159] ~ c(1:53), col = colIn[4], type = "b")


plot(state1Avg[54:106] ~ c(1:53), ylim = c(0, .6), type = "b",
     col = colIn[1], las = 1,
     xlab = "Week of year", ylab = "Pr[Animal in State]",
     main = "4 to 6.9 years")
points(state2Avg[54:106] ~ c(1:53), col = colIn[2], type = "b")
points(state3Avg[54:106] ~ c(1:53), col = colIn[3], type = "b")
points(state4Avg[54:106] ~ c(1:53), col = colIn[4], type = "b")


# by ewe repro status
fullTrajs$weekRepro <- paste(fullTrajs$ageCat, "_", formatC(fullTrajs$weekInYr, width = 2), sep = "")
state1Avg <- tapply(fullTrajs$state1,
                    fullTrajs$weekAge,
                    mean)
state2Avg <- tapply(fullTrajs$state2,
                    fullTrajs$weekAge,
                    mean)
state3Avg <- tapply(fullTrajs$state3,
                    fullTrajs$weekAge,
                    mean)
state4Avg <- tapply(fullTrajs$state4,
                    fullTrajs$weekAge,
                    mean)

alphaIn <- .5
colIn <- c(rgb(0, 0, 0, alpha = alphaIn),
           rgb(1, 0, 0, alpha = alphaIn),
           rgb(0, 1, 0, alpha = alphaIn),
           rgb(0, 0, 1, alpha = alphaIn),
           rgb(0, 1, 1, alpha = 1))
colIn <- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99')




par(mfrow = c(1, 3))
plot(state1Avg[160:212] ~ c(1:53), ylim = c(0, .6), type = "b",
     col = colIn[1], las = 1,
     xlab = "Day in year", ylab = "Pr[Animal in State]",
     main = "Over 6.9 years", xaxt = "n")
points(state2Avg[160:212] ~ c(1:53), col = colIn[2], type = "b")
points(state3Avg[160:212] ~ c(1:53), col = colIn[3], type = "b")
points(state4Avg[160:212] ~ c(1:53), col = colIn[4], type = "b")
axis(side = 1, at = c((31 + 28 + 31)/7,
                      (31 + 28 + 31 + 30 + 31 + 30)/7,
                      (31 + 28 + 31 + 30 + 31 + 30 + 31 + 31 + 30)/7),
     labels = c("Apr 01", "Jul 01", "Oct 01"))
abline(v = (31 + 28 + 31)/7, lwd = .5, 
       col = "grey40", lty = 2)
abline(v = (31 + 28 + 31 + 30 + 31 + 30)/7, lwd = .5, 
       col = "grey40", lty = 2)
abline(v = (31 + 28 + 31 + 30 + 31 + 30 + 31 + 31 + 30)/7, lwd = .5, 
       col = "grey40", lty = 2)

plot(state1Avg[107:159] ~ c(1:53), ylim = c(0, .6), type = "b",
     col = colIn[1], las = 1, xaxt = "n",
     xlab = "Day in year", ylab = "Pr[Animal in State]",
     main = "Under 4 years old")
points(state2Avg[107:159] ~ c(1:53), col = colIn[2], type = "b")
points(state3Avg[107:159] ~ c(1:53), col = colIn[3], type = "b")
points(state4Avg[107:159] ~ c(1:53), col = colIn[4], type = "b")
axis(side = 1, at = c((31 + 28 + 31)/7,
                      (31 + 28 + 31 + 30 + 31 + 30)/7,
                      (31 + 28 + 31 + 30 + 31 + 30 + 31 + 31 + 30)/7),
     labels = c("Apr 01", "Jul 01", "Oct 01"))
abline(v = (31 + 28 + 31)/7, lwd = .5, 
       col = "grey40", lty = 2)
abline(v = (31 + 28 + 31 + 30 + 31 + 30)/7, lwd = .5, 
       col = "grey40", lty = 2)
abline(v = (31 + 28 + 31 + 30 + 31 + 30 + 31 + 31 + 30)/7, lwd = .5, 
       col = "grey40", lty = 2)

plot(state1Avg[54:106] ~ c(1:53), ylim = c(0, .6), type = "b",
     col = colIn[1], las = 1, xaxt = "n", 
     xlab = "Day in year", ylab = "Pr[Animal in State]",
     main = "4 to 6.9 years")
points(state2Avg[54:106] ~ c(1:53), col = colIn[2], type = "b")
points(state3Avg[54:106] ~ c(1:53), col = colIn[3], type = "b")
points(state4Avg[54:106] ~ c(1:53), col = colIn[4], type = "b")
axis(side = 1, at = c((31 + 28 + 31)/7,
                      (31 + 28 + 31 + 30 + 31 + 30)/7,
                      (31 + 28 + 31 + 30 + 31 + 30 + 31 + 31 + 30)/7),
     labels = c("Apr 01", "Jul 01", "Oct 01"))
abline(v = (31 + 28 + 31)/7, lwd = .5, 
       col = "grey40", lty = 2)
abline(v = (31 + 28 + 31 + 30 + 31 + 30)/7, lwd = .5, 
       col = "grey40", lty = 2)
abline(v = (31 + 28 + 31 + 30 + 31 + 30 + 31 + 31 + 30)/7, lwd = .5, 
       col = "grey40", lty = 2)




modDirectedAvg <- tapply(fullTrajs$ModerateDirected,
                         fullTrajs$weekInYr,
                         mean)
staticAvg <- tapply(fullTrajs$StaticLocal,
                    fullTrajs$weekInYr,
                    mean)
modLocalAvg <- tapply(fullTrajs$ModerateLocal,
                      fullTrajs$weekInYr,
                      mean)
lowMovtLocalAvg <- tapply(fullTrajs$LowMovementStatic,
                          fullTrajs$weekInYr,
                          mean)

par(mfrow = c(1, 1))
plot(staticAvg ~ c(0:52), ylim = c(0, .5), type = "b", col = colIn[1])
points(lowMovtLocalAvg ~ c(0:52), col = colIn[2], type = "b")
points(modLocalAvg ~ c(0:52), col = colIn[3], type = "b")
points(modDirectedAvg ~ c(0:52), col = colIn[4], type = "b")

staticAgeAvg <- tapply(fullTrajs$StaticLocal,
                       fullTrajs$weekAge,
                       mean)
lowMoveLocalAgeAvg <- tapply(fullTrajs$LowMovementStatic,
                             fullTrajs$weekAge,
                             mean)
modDirAgeAvg <- tapply(fullTrajs$ModerateDirected,
                       fullTrajs$weekAge,
                       mean)
modLocalAgeAvg <- tapply(fullTrajs$ModerateLocal,
                         fullTrajs$weekAge,
                         mean)

par(mfrow = c(1, 3))
plot(staticAgeAvg[160:212] ~ c(1:53), ylim = c(0, .6), type = "b",
     col = colIn[1], las = 1, xaxt = "n", 
     xlab = "Day in year", ylab = "Pr[Animal in State]",
     main = "Under 4 years")
points(lowMoveLocalAgeAvg[160:212] ~ c(1:53), col = colIn[2], type = "b")
points(modLocalAgeAvg[160:212] ~ c(1:53), col = colIn[3], type = "b")
points(modDirAgeAvg[160:212] ~ c(1:53), 
       col = colIn[4], type = "b", pch = 16, lwd = 3)
axis(side = 1, at = c((31 + 28 + 31)/7,
                      (31 + 28 + 31 + 30 + 31 + 30)/7,
                      (31 + 28 + 31 + 30 + 31 + 30 + 31 + 31 + 30)/7),
     labels = c("Apr 01", "Jul 01", "Oct 01"))
abline(v = (31 + 28 + 31)/7, lwd = .5, 
       col = "grey40", lty = 2)
abline(v = (31 + 28 + 31 + 30 + 31 + 30)/7, lwd = .5, 
       col = "grey40", lty = 2)
abline(v = (31 + 28 + 31 + 30 + 31 + 30 + 31 + 31 + 30)/7, lwd = .5, 
       col = "grey40", lty = 2)

plot(staticAgeAvg[107:159] ~ c(1:53), ylim = c(0, .6), type = "b",
     col = colIn[1], las = 1, xaxt = "n",
     xlab = "Day in year", ylab = "Pr[Animal in State]",
     main = "4.0 - 6.9 years")
points(lowMoveLocalAgeAvg[107:159] ~ c(1:53), col = colIn[2], type = "b")
points(modLocalAgeAvg[107:159] ~ c(1:53), col = colIn[3], type = "b")
points(modDirAgeAvg[107:159] ~ c(1:53), 
       col = colIn[4], type = "b", pch = 16, lwd = 3)
axis(side = 1, at = c((31 + 28 + 31)/7,
                      (31 + 28 + 31 + 30 + 31 + 30)/7,
                      (31 + 28 + 31 + 30 + 31 + 30 + 31 + 31 + 30)/7),
     labels = c("Apr 01", "Jul 01", "Oct 01"))
abline(v = (31 + 28 + 31)/7, lwd = .5, 
       col = "grey40", lty = 2)
abline(v = (31 + 28 + 31 + 30 + 31 + 30)/7, lwd = .5, 
       col = "grey40", lty = 2)
abline(v = (31 + 28 + 31 + 30 + 31 + 30 + 31 + 31 + 30)/7, lwd = .5, 
       col = "grey40", lty = 2)

plot(staticAgeAvg[54:106] ~ c(1:53), ylim = c(0, .6), type = "b",
     col = colIn[1], las = 1, xaxt = "n", 
     xlab = "Day in year", ylab = "Pr[Animal in State]",
     main = "Over 6.9 years")
points(lowMoveLocalAgeAvg[54:106] ~ c(1:53), col = colIn[2], type = "b")
points(modLocalAgeAvg[54:106] ~ c(1:53), col = colIn[3], type = "b")
points(modDirAgeAvg[54:106] ~ c(1:53), 
       col = colIn[4], type = "b", lwd = 3, pch = 16)
axis(side = 1, at = c((31 + 28 + 31)/7,
                      (31 + 28 + 31 + 30 + 31 + 30)/7,
                      (31 + 28 + 31 + 30 + 31 + 30 + 31 + 31 + 30)/7),
     labels = c("Apr 01", "Jul 01", "Oct 01"))
abline(v = (31 + 28 + 31)/7, lwd = .5, 
       col = "grey40", lty = 2)
abline(v = (31 + 28 + 31 + 30 + 31 + 30)/7, lwd = .5, 
       col = "grey40", lty = 2)
abline(v = (31 + 28 + 31 + 30 + 31 + 30 + 31 + 31 + 30)/7, lwd = .5, 
       col = "grey40", lty = 2)

