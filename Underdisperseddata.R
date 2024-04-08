# Clear the workspace
rm(list = ls())

# Set seed for reproducibility
set.seed(123)

# Load required libraries
library(hnp)
library(tidyverse)
library(gamlss)
library(GLDEX)
library(pscl)

# Define a function to suppress output
quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
} 

# Define a function to calculate distance
my_distance <- function(hnp_obj, f = function(w) 1, g = function(x, alpha = 2, delta = 0.2) 1) {
  r <- hnp_obj$residuals
  m <- hnp_obj$median
  u <- hnp_obj$upper
  l <- hnp_obj$lower
  Indicator <- as.numeric(r > u | r < l)
  w <- u - l
  
  b <- ifelse(r > u, r - u, ifelse(r < l, l - r, 1))
  distance <- sum((abs(r - m) / f(w)) * g(b) ^ Indicator)
  distance_w <- sum(((r - m) ^ 2 / f(w)))
  my_list <- list(distance, f(w), g(b), sum(abs(r - m)), w, b, r, m, distance_w)
  
  return(my_list)
}

# Define the main simulation function
func_sim <- function(n, z, beta_coef = c(3, 2.5), theta, f = function(w) 1, g = function(x, alpha = 2, delta = 0.2) 1) {
  x <- rnorm(n)
  mu <- exp(beta_coef[1] + beta_coef[2] * x)
  y <- t(replicate(z, 0.5 * rpois(n, mu)))
  
  # Initialize lists to store results
  fit_poisson <- fit_quasi <- my_hnp_poisson <- my_hnp_quasi <- list()
  distance_poisson <- distance_quasi <- numeric(z)
  distance_w_poisson <- distance_w_quasi <- numeric(z)
  b_poisson <- w_poisson <- b_quasi <- w_quasi <- numeric(z)
  base_value_poisson <- base_value_quasi <- numeric(z)
  w_value_poisson <- w_value_quasi <- numeric(z)
  b_value_poisson <- b_value_quasi <- numeric(z)
  r_poisson <- r_quasi <- numeric(z)
  m_poisson <- m_quasi <- numeric(z)
  aic_poisson <- bic_poisson <- numeric(z)
  
  # Loop through simulations
  for (i in 1:z) {
    # Fit Poisson and Quasi-Poisson models
    fit_poisson[[i]] <- glm(y[i,] ~ x, family = "poisson")
    fit_quasi[[i]] <- glm(y[i,] ~ x, family = "quasipoisson")
    
    # Perform heteroscedasticity-robust hypothesis testing
    my_hnp_poisson[[i]] <- quiet(hnp(fit_poisson[[i]], plot = FALSE, resid.type = "pearson"))
    my_hnp_quasi[[i]] <- quiet(hnp(fit_quasi[[i]], plot = FALSE, resid.type = "pearson"))
    
    # Calculate distances and related values
    distance_poisson[i] <- my_distance(my_hnp_poisson[[i]], f = f, g = g)[[1]]
    distance_quasi[i] <- my_distance(my_hnp_quasi[[i]], f = f, g = g)[[1]]
    distance_w_poisson[i] <- my_distance(my_hnp_poisson[[i]], f = f, g = g)[[9]]
    distance_w_quasi[i] <- my_distance(my_hnp_quasi[[i]], f = f, g = g)[[9]]
    b_poisson[i] <- my_distance(my_hnp_poisson[[i]], f = f, g = g)[[3]]
    w_poisson[i] <- my_distance(my_hnp_poisson[[i]], f = f, g = g)[[2]]
    b_quasi[i] <- my_distance(my_hnp_quasi[[i]], f = f, g = g)[[3]]
    w_quasi[i] <- my_distance(my_hnp_quasi[[i]], f = f, g = g)[[2]]
    base_value_poisson[i] <- my_distance(my_hnp_poisson[[i]], f = f, g = g)[[4]]
    base_value_quasi[i] <- my_distance(my_hnp_quasi[[i]], f = f, g = g)[[4]]
    w_value_poisson[i] <- my_distance(my_hnp_poisson[[i]], f = f, g = g)[[5]]
    w_value_quasi[i] <- my_distance(my_hnp_quasi[[i]], f = f, g = g)[[5]]
    b_value_poisson[i] <- my_distance(my_hnp_poisson[[i]], f = f, g = g)[[6]]
    b_value_quasi[i] <- my_distance(my_hnp_quasi[[i]], f = f, g = g)[[6]]
    r_poisson[i] <- my_distance(my_hnp_poisson[[i]], f = f, g = g)[[7]]
    r_quasi[i] <- my_distance(my_hnp_quasi[[i]], f = f, g = g)[[7]]
    m_poisson[i] <- my_distance(my_hnp_poisson[[i]], f = f, g = g)[[8]]
    m_quasi[i] <- my_distance(my_hnp_quasi[[i]], f = f, g = g)[[8]]
    aic_poisson[i] <- AIC(fit_poisson[[i]])
    bic_poisson[i] <- BIC(fit_poisson[[i]])
  }
  
  # Store results in a list
  my_return <- list("model_poisson" = fit_poisson, "model_quasi" = fit_quasi,
                    "hnp_poisson" = my_hnp_poisson, "hnp_quasi" = my_hnp_quasi,
                    "dist_poisson" = distance_poisson, "dist_quasi" = distance_quasi,
                    "dist_w_poisson" = distance_w_poisson, "dist_w_quasi" = distance_w_quasi,
                    "b_poisson" = b_poisson, "w_poisson" = w_poisson, "b_quasi" = b_quasi, "w_quasi" = w_quasi,
                    "base_value_poisson" = base_value_poisson, "base_value_quasi" = base_value_quasi,
                    "w_value_poisson" = w_value_poisson, "w_value_quasi" = w_value_quasi,
                    "b_value_poisson" = b_value_poisson, "b_value_quasi" = b_value_quasi,
                    "r_poisson" = r_poisson, "r_quasi" = r_quasi,
                    "m_poisson" = m_poisson, "m_quasi" = m_quasi,
                    "aic_poisson" = aic_poisson, "bic_poisson" = bic_poisson)
  
  return(invisible(my_return))
}

# Run the simulation
my_simulation1 <- func_sim(n = 20, z = 1000)

# Save the results
save(my_simulation1, file = "sim_upois1_new.Rdata")
