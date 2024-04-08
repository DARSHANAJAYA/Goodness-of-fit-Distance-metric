
# Clearing the current environment
rm(list = ls())

# Setting seed for reproducibility
set.seed(123)

# Loading required libraries
library(hnp)      # Provides functions for Heteroscedastic Nonlinear Model-based Regression
library(tidyverse)  # Collection of packages for data manipulation and visualization
library(gamlss)    # Generalized Additive Models for Location Scale and Shape
library(GLDEX)     # Provides functions to fit generalized lambda distributions
library(pscl)      # Provides functions for modeling count data regression models

# Function to suppress output to the console
quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
} 

# Function to calculate distance metric for heteroscedasticity
my_distance <- function(hnp_obj, f = function(w) 1, g = function(x, alpha = 2, delta = 0.2) 1) {
  r <- hnp_obj$residuals  # Residuals
  m <- hnp_obj$median      # Median
  u <- hnp_obj$upper       # Upper bound
  l <- hnp_obj$lower       # Lower bound
  
  # Indicator function for residuals outside of bounds
  Indicator <- as.numeric(r > u | r < l)
  
  # Weight calculation
  w <- u - l
  
  # Calculation of distance
  b <- ifelse(r > u, r - u, ifelse(r < l, l - r, 1))
  distance <- sum(((r - m)^2 / f(w)) * g(b)^Indicator)
  distance_w <- sum(((r - m)^2 / f(w)))
  
  my_list <- list(distance, f(w), g(b), sum((r - m)^2), w, b, r, m, distance_w)
  
  return(my_list)
}

# Function to simulate data and fit models
func_sim <- function(n, z, beta_coef = c(3, 2.5), theta, f = function(w) 1, g = function(x, alpha = 2, delta = 0.2) 1) {
  x <- rnorm(n)  # Generating input data
  
  # Calculating mean using linear model
  mu <- exp(beta_coef[1] + beta_coef[2] * x)
  
  # Generating count data
  y <- t(replicate(z, 0.5 * rpois(n, mu)))
  
  # Fitting various models and calculating distances
  # Looping over simulations
  for (i in 1:z) {
    fit_nb_lin[[i]] <- try(quiet(gamlss(y[i, ] ~ x, family = NBII)), silent = TRUE)
    fit_nb_quad[[i]] <- try(glm.nb(y[i, ] ~ x), silent = TRUE)
    
    if (class(fit_nb_lin[[i]]) == "try-error" | class(fit_nb_quad[[i]]) == "try-error") {
      next
    }
    
    if (fit_nb_lin[[i]]["converged"] == FALSE) {
      next
    }
    
    if (fit_nb_lin[[i]]$sigma.fv[1] > 1000) {
      next
    }
    
    fit_poisson[[i]] <- glm(y[i, ] ~ x, family = "poisson")
    fit_quasi[[i]] <- glm(y[i, ] ~ x, family = "quasipoisson")
    
    my_hnp_nb_quad[[i]] <- quiet(hnp(fit_nb_quad[[i]], plot = FALSE, resid.type = "pearson"))
    my_hnp_poisson[[i]] <- quiet(hnp(fit_poisson[[i]], plot = FALSE, resid.type = "pearson"))
    my_hnp_quasi[[i]] <- quiet(hnp(fit_quasi[[i]], plot = FALSE, resid.type = "pearson"))
    
    my_hnp_nb_lin[[i]] <- quiet(hnp(fit_nb_lin[[i]], newclass = TRUE, resid.type = "pearson", plot = FALSE))
    
    # Calculating distances
    distance_nb_quad[i] <- my_distance(my_hnp_nb_quad[[i]], f = f, g = g)[[1]]
    distance_poisson[i] <- my_distance(my_hnp_poisson[[i]], f = f, g = g)[[1]]
    distance_quasi[i] <- my_distance(my_hnp_quasi[[i]], f = f, g = g)[[1]]
    distance_nb_lin[i] <- my_distance(my_hnp_nb_lin[[i]], f = f, g = g)[[1]]
    
    distance_w_nb_quad[i] <- my_distance(my_hnp_nb_quad[[i]], f = f, g = g)[[9]]
    distance_w_poisson[i] <- my_distance(my_hnp_poisson[[i]], f = f, g = g)[[9]]
    distance_w_quasi[i] <- my_distance(my_hnp_quasi[[i]], f = f, g = g)[[9]]
    distance_w_nb_lin[i] <- my_distance(my_hnp_nb_lin[[i]], f = f, g = g)[[9]]
    
    # Other metrics
    # ...
  }
  
  # Constructing return object
  my_return <- list(
    # Models
    "model_nb_quad" = fit_nb_quad,
    "model_poisson" = fit_poisson,
    "model_quasi" = fit_quasi,
    "model_nb_lin" = fit_nb_lin,
    
    # HNP objects
    "hnp_nb_quad" = my_hnp_nb_quad,
    "hnp_nb_lin" = my_hnp_nb_lin,
    "hnp_poisson" = my_hnp_poisson,
    "hnp_quasi" = my_hnp_quasi,
    
    # Distances
    "dist_nb_lin" = distance_nb_lin,
    "dist_nb_quad" = distance_nb_quad,
    "dist_poisson" = distance_poisson,
    "dist_quasi" = distance_quasi,
    
    # Weighted distances
    "dist_w_nb_lin" = distance_w_nb_lin,
    "dist_w_nb_quad" = distance_w_nb_quad,
    "dist_w_poisson" = distance_w_poisson,
    "dist_w_quasi" = distance_w_quasi,
    
    # Other metrics
    # ...
  )
  
  return(invisible(my_return))
}

# Generating simulated data and saving results
set.seed(123)
my_simulation1 <- func_sim(n = 20, z = 1000)
save(my_simulation1, file = "sim_pois1_new.Rdata")