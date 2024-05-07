# Clear the environment
rm(list = ls())

# Set seed for reproducibility
set.seed(123)

# Load required libraries
library(hnp)   # Provides functions for heteroscedastic non-parametric tests
library(tidyverse)  # For data manipulation and visualization
library(gamlss)   # Generalized Additive Models for Location, Scale and Shape
library(GLDEX)   # Generalized Lambda Distributions and Extreme Value Distributions

# Define a function to suppress output
quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
} 

# Function to calculate distance
my_distance <- function(hnp_obj, f = function(w) 1, g = function(x, alpha = 2, delta = 0.2) 1) {
  r <- hnp_obj$residuals
  m <- hnp_obj$median
  u <- hnp_obj$upper
  l <- hnp_obj$lower
  Indicator <- as.numeric(r > u | r < l)  # Indicator function
  
  w <- u - l  # Width of confidence interval
  
  # Calculating the distance based on provided formulas
  distance <- sum(((r - m)^2 / f(w)) * g(ifelse(r > u, r - u, ifelse(r < l, l - r, 1)))^Indicator)
  
  # Return a list containing relevant values
  my_list <- list(distance, f(w), g(ifelse(r > u, r - u, ifelse(r < l, l - r, 1))), sum((r - m)^2), w, ifelse(r > u, r - u, ifelse(r < l, l - r, 1)), r, m)
  
  return(my_list)
}

# Function to simulate data
func_sim <- function(n, z, beta_coef = c(3, 2), sigma, f = function(w) 1, g = function(x, alpha = 2, delta = 0.2) 1) {
  x <- rnorm(n)
  mu <- exp(beta_coef[1] + beta_coef[2] * x)  # Generate mean values
  
  # Simulate data using Negative Binomial II distribution
  y <- t(replicate(z, rNBII(n, mu = mu, sigma)))
  
  # Initialize empty lists for storing results
  fit_nb_quad <- fit_nb_lin <- fit_poisson <- fit_quasi <- my_hnp_nb_quad <- my_hnp_nb_lin <- my_hnp_poisson <- my_hnp_quasi <- b_nb_quad <- b_nb_lin <- b_poisson <- b_quasi <- w_nb_quad <- w_nb_lin <- w_poisson <- w_quasi <- aic_neg_lin <- aic_neg_quad <- aic_poisson <- bic_neg_lin <- bic_neg_quad <- bic_poisson <- list()
  base_value_poisson <- base_value_negbinlin <- base_value_negbinquad <- base_value_quasi <- w_value_poisson <- w_value_negbinlin <- w_value_negbinquad <- w_value_quasi <- b_value_poisson <- b_value_negbinlin <- b_value_negbinquad <- b_value_quasi <- list()
  m_poisson <- m_negbinquad <- m_negbinlin <- m_quasi <- list()
  r_poisson <- r_negbinquad <- r_negbinlin <- r_quasi <- list()
  distance_nb_lin <- distance_nb_quad <- distance_poisson <- distance_quasi <- numeric(z)
  
  # Define a function to calculate residuals
  dfun <- function(obj) {
    y <- obj$y
    mu_hat <- predict(fit_nb_lin[[i]], what = "mu", type = "response")
    sigma_hat <- predict(fit_nb_lin[[i]], what = "sigma", type = "response")
    v_mu <- (1 + sigma_hat) * mu_hat
    res <- (y - mu_hat) / sqrt(v_mu)
    return(res)
  }
  
  # Define a function to simulate response data
  sfun <- function(n, obj) {
    mu_hat <- predict(fit_nb_lin[[i]], what = "mu", type = "response")
    sigma_hat <- predict(fit_nb_lin[[i]], what = "sigma", type = "response")
    sim <- rNBII(n = length(mu_hat), mu = mu_hat, sigma = sigma_hat)
    return(sim)
  }
  
  # Define a function to fit models
  ffun <- function(resp) {
    try_fit <- try(gamlss(resp ~ x, family = NBII), silent = TRUE)
    n_try <- 1
    while(class(try_fit) == "try-error" & n_try < 20) {
      new_response <- sfun(1, object)
      try_fit <- try(gamlss(new_response ~ x, family = NBII), silent = TRUE)
      n_try <- n_try + 1
    }
    return(try_fit)
  }
  
  # Loop over simulations
  for (i in 1:z) {
    fit_nb_lin[[i]] <- try(quiet(gamlss(y[i,] ~ x, family = NBII)), silent = TRUE)
    fit_nb_quad[[i]] <- try(glm.nb(y[i,] ~ x), silent = TRUE)
    
    # Skip to the next iteration if fitting fails or conditions are not met
    if (class(fit_nb_lin[[i]]) == "try-error" | class(fit_nb_quad[[i]]) == "try-error" | fit_nb_lin[[i]]["converged"] == FALSE | fit_nb_lin[[i]]$sigma.fv[1] > 1000) {
      next
    }
    
    # Fit Poisson and Quasi-Poisson models
    fit_poisson[[i]] <- glm(y[i,] ~ x, family = "poisson")
    fit_quasi[[i]] <- glm(y[i,] ~ x, family = "quasipoisson")
    
    # Compute HNP (Heteroscedastic Non-parametric test) for each model
    my_hnp_nb_quad[[i]] <- quiet(hnp(fit_nb_quad[[i]], plot = FALSE, resid.type = "pearson"))
    my_hnp_poisson[[i]] <- quiet(hnp(fit_poisson[[i]], plot = FALSE, resid.type = "pearson"))
    my_hnp_quasi[[i]] <- quiet(hnp(fit_quasi[[i]], plot = FALSE, resid.type = "pearson"))
    my_hnp_nb_lin[[i]] <- quiet(hnp(fit_nb_lin[[i]], newclass = TRUE, resid.type = "pearson", diagfun = dfun, simfun = sfun, fitfun = ffun, plot = FALSE))
    
    # Calculate distances using custom function
    distance_nb_quad[i] <- my_distance(my_hnp_nb_quad[[i]], f = f, g = g)[[1]]
    distance_poisson[i] <- my_distance(my_hnp_poisson[[i]], f = f, g = g)[[1]]
    distance_quasi[i] <- my_distance(my_hnp_quasi[[i]], f = f, g = g)[[1]]
    distance_nb_lin[i] <- my_distance(my_hnp_nb_lin[[i]], f = f, g = g)[[1]]
    
    # Store additional parameters
    b_nb_quad[i] <- my_distance(my_hnp_nb_quad[[i]], f = f, g = g)[[3]]
    w_nb_quad[i] <- my_distance(my_hnp_nb_quad[[i]], f = f, g = g)[[2]]
    b_poisson[i] <- my_distance(my_hnp_poisson[[i]], f = f, g = g)[[3]]
    w_poisson[i] <- my_distance(my_hnp_poisson[[i]], f = f, g = g)[[2]]
    b_quasi[i] <- my_distance(my_hnp_quasi[[i]], f = f, g = g)[[3]]
    w_quasi[i] <- my_distance(my_hnp_quasi[[i]], f = f, g = g)[[2]]
    b_nb_lin[i] <- my_distance(my_hnp_nb_lin[[i]], f = f, g = g)[[3]]
    w_nb_lin[i] <- my_distance(my_hnp_nb_lin[[i]], f = f, g = g)[[2]]
    base_value_poisson[i] <- my_distance(my_hnp_poisson[[i]], f = f, g = g)[[4]]
    base_value_negbinlin[i] <- my_distance(my_hnp_nb_lin[[i]], f = f, g = g)[[4]]
    base_value_negbinquad[i] <- my_distance(my_hnp_nb_quad[[i]], f = f, g = g)[[4]] 
    base_value_quasi[i] <- my_distance(my_hnp_quasi[[i]], f = f, g = g)[[4]]
    w_value_poisson[i] <- my_distance(my_hnp_poisson[[i]], f = f, g = g)[[5]]
    w_value_negbinlin[i] <- my_distance(my_hnp_nb_lin[[i]], f = f, g = g)[[5]]
    w_value_negbinquad[i] <- my_distance(my_hnp_nb_quad[[i]], f = f, g = g)[[5]] 
    w_value_quasi[i] <- my_distance(my_hnp_quasi[[i]], f = f, g = g)[[5]]
    b_value_poisson[i] <- my_distance(my_hnp_poisson[[i]], f = f, g = g)[[6]]
    b_value_negbinlin[i] <- my_distance(my_hnp_nb_lin[[i]], f = f, g = g)[[6]]
    b_value_negbinquad[i] <- my_distance(my_hnp_nb_quad[[i]], f = f, g = g)[[6]] 
    b_value_quasi[i] <- my_distance(my_hnp_quasi[[i]], f = f, g = g)[[6]]
    r_poisson[i] <- my_distance(my_hnp_poisson[[i]], f = f, g = g)[[7]]
    r_negbinlin[i] <- my_distance(my_hnp_nb_lin[[i]], f = f, g = g)[[7]]
    r_negbinquad[i] <- my_distance(my_hnp_nb_quad[[i]], f = f, g = g)[[7]] 
    r_quasi[i] <- my_distance(my_hnp_quasi[[i]], f = f, g = g)[[7]]
    m_poisson[i] <- my_distance(my_hnp_poisson[[i]], f = f, g = g)[[8]]
    m_negbinlin[i] <- my_distance(my_hnp_nb_lin[[i]], f = f, g = g)[[8]]
    m_negbinquad[i] <- my_distance(my_hnp_nb_quad[[i]], f = f, g = g)[[8]] 
    m_quasi[i] <- my_distance(my_hnp_quasi[[i]], f = f, g = g)[[8]]
    aic_neg_lin[i] <- AIC(fit_nb_lin[[i]])
    aic_neg_quad[i] <- AIC(fit_nb_quad[[i]])
    aic_poisson[i] <- AIC(fit_poisson[[i]])
    bic_neg_lin[i] <- BIC(fit_nb_lin[[i]])
    bic_neg_quad[i] <- BIC(fit_nb_quad[[i]])
    bic_poisson[i] <- BIC(fit_poisson[[i]])
  }
  
  # Compute summary statistics for distances
  summary_distance <- data.frame("distance_nb" = c(summary(distance_nb_quad), sd(distance_nb_quad)),
                                 "distance_poisson" = c(summary(distance_poisson), sd(distance_poisson)),
                                 "distance_quasi" = c(summary(distance_quasi), sd(distance_quasi)),
                                 "distance_nb_lin" = c(summary(distance_nb_lin), sd(distance_nb_lin)))
  row.names(summary_distance)[7] <- "std_deviation"
  
  print(summary_distance)
  
  # Store results in a list
  my_return <- list("model_nb_quad" = fit_nb_quad, "model_poisson" = fit_poisson, "model_quasi" = fit_quasi, "model_nb_lin" = fit_nb_lin,
                    "hnp_nb_quad" = my_hnp_nb_quad, "hnp_nb_lin" = my_hnp_nb_lin,
                    "hnp_poisson" = my_hnp_poisson, "hnp_quasi" = my_hnp_quasi,
                    "dist_nb_lin" = distance_nb_lin, "dist_nb_quad" = distance_nb_quad,
                    "dist_poisson" = distance_poisson, "dist_quasi" = distance_quasi,
                    "b_nb_quad" = b_nb_quad, "w_nb_quad" = w_nb_quad, "b_poisson" = b_poisson, "w_poisson" = w_poisson,
                    "b_quasi" = b_quasi, "w_quasi" = w_quasi, "b_nb_lin" = b_nb_lin, "w_nb_lin" = w_nb_lin,
                    "aic_neg_lin" = aic_neg_lin, "aic_neg_quad" = aic_neg_quad, "aic_poisson" = aic_poisson,
                    "bic_neg_lin" = bic_neg_lin, "bic_neg_quad" = bic_neg_quad, "bic_poisson" = bic_poisson, 
                    "base_value_poisson" = base_value_poisson, "base_value_negbinlin" = base_value_negbinlin,
                    "base_value_negbinquad" = base_value_negbinquad, "base_value_quasi" = base_value_quasi,
                    "w_value_poisson" = w_value_poisson, "w_value_negbinlin" = w_value_negbinlin,
                    "w_value_negbinquad" = w_value_negbinquad, "w_value_quasi" = w_value_quasi,
                    "b_value_poisson" = b_value_poisson, "b_value_negbinlin" = b_value_negbinlin,
                    "b_value_negbinquad" = b_value_negbinquad, "b_value_quasi" = b_value_quasi,
                    "r_poisson" = r_poisson, "r_negbinquad" = r_negbinquad, "r_negbinlin" = r_negbinlin,
                    "r_quasi" = r_quasi, "m_poisson" = m_poisson, "m_negbinquad" = m_negbinquad,
                    "m_negbinlin" = m_negbinlin, "m_quasi" = m_quasi, "summary" = summary_distance)
  
  return(invisible(my_return))
}

# Set seed and run simulation
set.seed(123)
my_simulation1 <- func_sim(n = 20, z = 1000, sigma = 7)

# Save the results
save(my_simulation1, file = "sim_negbinhlin1_new.Rdata")
