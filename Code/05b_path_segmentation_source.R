
# I. gamma_pars ----
## This function is from Simona Picardi. It converts a shape/scale
## parameterization of the gamma to a mean/standard deviation parameterization.
## (This is essentially achieved with moment-matching). 
gamma_pars <- function(mean = NULL, sd = NULL,
                       rate = NULL, shape = NULL, scale = NULL) {
  # case 1: if mean and sd are provided, convert to shape and scale
  if(is.numeric(mean) && is.numeric(sd)) {
    
    rate <- mean/sd^2
    shape <- (mean/sd)^2
    
    pars <- c(rate = rate,
              shape = shape)
    
    # case 2: if rate and shape are given, convert to mean and sd
  } else if(is.numeric(rate) && is.numeric(shape)) {
    
    mean <- shape/rate
    sd <- sqrt(shape)/rate
    
    pars <- c(mean = mean,
              sd = sd)
    
    # case 3: if shape and scale are given, convert to mean and sd
  } else if(is.numeric(scale) && is.numeric(shape)) {
    
    rate <- 1/scale
    
    mean <- shape/rate
    sd <- sqrt(shape)/rate
    
    pars <- c(mean = mean,
              sd = sd)
    
  }
  return(pars)
}
