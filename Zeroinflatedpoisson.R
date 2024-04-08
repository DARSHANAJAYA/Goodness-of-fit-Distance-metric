# Clear the workspace
rm(list = ls())

# Set seed for reproducibility
set.seed(5678)

# Load required libraries
library(hnp)
library(tidyverse)
library(gamlss)
library(GLDEX)
library(pscl)
library(data.table)

# Define a function to suppress output
quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
} 

# Define a custom distance function
my_distance <- function(hnp_obj, f = function(w) 1, g = function(x, alpha = 2, delta = 0.2) 1) {
  r <- hnp_obj$residuals
  m <- hnp_obj$median
  u <- hnp_obj$upper
  l <- hnp_obj$lower
  Indicator <- as.numeric(r > u | r < l)
  w <- u - l
  
  b <- ifelse(r > u, r - u, ifelse(r < l, l - r, 1))
  distance <- sum(((r - m)^2 / f(w)) * g(b)^Indicator)
  distance_w <- sum(((r - m)^2 / f(w)))
  my_list <- list(distance, f(w), g(b), sum((r - m)^2), w, b, r, m, distance_w)
  
  return(my_list)
}

# Define the main simulation function
func_sim <- function(n, z, beta_coef = c(3, 2.5), theta, f = function(w) 1, g = function(x, alpha = 2, delta = 0.2) 1) {
  x <- rnorm(n)
  mu <- exp(beta_coef[1] + beta_coef[2] * x)
  
  # Simulate data using rZINBI
  y <- t(replicate(z, rZINBI(n, mu, sigma = 7, nu = 0.6))) # sigma = 7, overdispersed; nu = probability of a zero obs
  
  # Initialize lists to store results
  fit_nb_quad <- fit_nb_lin <- fit_poisson <- fit_quasi <- fit_zinbh <- fit_zinp <- my_hnp_nb_quad <- my_hnp_nb_lin <- my_hnp_poisson <- my_hnp_quasi <- my_hnp_zinbh <- my_hnp_zinp <- b_nb_quad <- b_nb_lin <- b_poisson <- b_quasi <- b_zinbh <- b_zinp <- w_nb_quad <- w_nb_lin <- w_poisson <- w_quasi <- w_zinbh <- w_zinp <- aic_neg_lin <- aic_neg_quad <- aic_poisson <- aic_zinbh <- aic_zinp <- bic_neg_lin <- bic_neg_quad <- bic_poisson <- bic_zinbh <- bic_zinp <- list()
  base_value_poisson <- base_value_negbinlin <- base_value_negbinquad <- base_value_quasi <- base_value_zinbh <- base_value_zinp <- w_value_poisson <- w_value_negbinlin <- w_value_negbinquad <- w_value_quasi <- w_value_zinbh <- w_value_zinp <- b_value_poisson <- b_value_negbinlin <- b_value_negbinquad <- b_value_quasi <- b_value_zinbh <- b_value_zinp <- list()
  m_poisson <- m_negbinquad <- m_negbinlin <- m_quasi <- m_zinbh <- m_zinp <- list()
  r_poisson <- r_negbinquad <- r_negbinlin <- r_quasi <- r_zinbh <- r_zinp <- list()
  distance_nb_lin <- distance_nb_quad <- distance_poisson <- distance_quasi <- distance_zinbh <- distance_zinp <- numeric(z)
  distance_w_nb_lin <- distance_w_nb_quad <- distance_w_poisson <- distance_w_quasi <- distance_w_zinbh <- distance_w_zinp <- numeric(z)
  nu_hat <- list()
  
  # Define a function to calculate residuals for NBII model
  dfun <- function(obj) {
    y <- obj$y
    mu_hat <- predict(fit_nb_lin[[i]], what = "mu", type = "response")
    sigma_hat <- predict(fit_nb_lin[[i]], what = "sigma", type = "response")
    v_mu <- (1 + sigma_hat) * mu_hat
    res <- (y - mu_hat) / sqrt(v_mu)
    return(res)
  }
  
  # Define a function to simulate data for NBII model
  sfun <- function(n, obj) {
    mu_hat <- predict(fit_nb_lin[[i]], what = "mu", type = "response")
    sigma_hat <- predict(fit_nb_lin[[i]], what = "sigma", type = "response")
    sim <- rNBII(n = length(mu_hat), mu = mu_hat, sigma = sigma_hat)
    return(sim)
  }
  
  # Define a function to fit NBII model
  ffun <- function(resp) {
    try_fit <- try(gamlss(resp ~ x, family = NBII), silent = TRUE)
    n_try <- 1
    while (class(try_fit) == "try-error" & n_try < 20) {
      new_response <- sfun(1, object)
      try_fit <- try(gamlss(new_response ~ x, family = NBII), silent = TRUE)
      n_try <- n_try + 1
    }
    return(try_fit)
  }
  
  # Similar functions for ZIP and ZINBI models...
  
  # Loop through simulations
  for (i in 1:z) {
    # Fit models
    fit_zinbh[[i]] <- try(quiet(zeroinfl(y[i,] ~ x | 1, dist = "negbin")), silent = TRUE)
    fit_zinp[[i]] <- try(quiet(zeroinfl(y[i,] ~ x | 1, dist = "poisson")), silent = TRUE)
    fit_nb_lin[[i]] <- try(quiet(gamlss(y[i,] ~ x, family = NBII)), silent = TRUE)
    fit_nb_quad[[i]] <- try(glm.nb(y[i,] ~ x), silent = TRUE)
    
    # Check for errors or convergence issues
    if (class(fit_nb_lin[[i]]) == "try-error" | class(fit_nb_quad[[i]]) == "try-error") {
      next
    }
    if (!fit_nb_lin[[i]]$converged) {
      next
    }
    if (fit_nb_lin[[i]]$sigma.fv[1] > 1000) {
      next
    }
    
    # Fit Poisson and Quasi-Poisson models
    fit_poisson[[i]] <- glm(y[i,] ~ x, family = "poisson")
    fit_quasi[[i]] <- glm(y[i,] ~ x, family = "quasipoisson")
    
    # Fit zero-inflated negative binomial and Poisson models
    my_hnp_nb_quad[[i]] <- quiet(hnp(fit_nb_quad[[i]], plot = FALSE, resid.type = "pearson"))
    my_hnp_poisson[[i]] <- quiet(hnp(fit_poisson[[i]], plot = FALSE, resid.type = "pearson"))
    my_hnp_quasi[[i]] <- quiet(hnp(fit_quasi[[i]], plot = FALSE, resid.type = "pearson"))
    my_hnp_zinbh[[i]] <- quiet(hnp(fit_zinbh[[i]], newclass = TRUE, diagfun = dfun_zinb, simfun = sfun_zinb, fitfun = ffun_zinb, plot = FALSE))
    my_hnp_zinp[[i]] <- quiet(hnp(fit_zinp[[i]], newclass = TRUE, diagfun = dfun_zip, simfun = sfun_zip, fitfun = ffun_zip, plot = FALSE))
    my_hnp_nb_lin[[i]] <- quiet(hnp(fit_nb_lin[[i]], newclass = TRUE, diagfun = dfun, simfun = sfun, fitfun = ffun, plot = FALSE))
    
    # Calculate distances
    distance_nb_quad[i] <- my_distance(my_hnp_nb_quad[[i]], f = f, g = g)[[1]]
    distance_poisson[i] <- my_distance(my_hnp_poisson[[i]], f = f, g = g)[[1]]
    distance_quasi[i] <- my_distance(my_hnp_quasi[[i]], f = f, g = g)[[1]]
    distance_nb_lin[i] <- my_distance(my_hnp_nb_lin[[i]], f = f, g = g)[[1]]
    distance_zinbh[i] <- my_distance(my_hnp_zinbh[[i]], f = f, g = g)[[1]]
    distance_zinp[i] <- my_distance(my_hnp_zinp[[i]], f = f, g = g)[[1]]
    
    distance_w_nb_quad[i] <- my_distance(my_hnp_nb_quad[[i]], f = f, g = g)[[9]]
    distance_w_poisson[i] <- my_distance(my_hnp_poisson[[i]], f = f, g = g)[[9]]
    distance_w_quasi[i] <- my_distance(my_hnp_quasi[[i]], f = f, g = g)[[9]]
    distance_w_nb_lin[i] <- my_distance(my_hnp_nb_lin[[i]], f = f, g = g)[[9]]
    distance_w_zinbh[i] <- my_distance(my_hnp_zinbh[[i]], f = f, g = g)[[9]]
    distance_w_zinp[i] <- my_distance(my_hnp_zinp[[i]], f = f, g = g)[[9]]
    
    # Extract other values of interest
    # ...
    
    cat("Simulation number:", i, "\n")
  }
  
  # Return results
  my_return <- list("model_nb_quad"= fit_nb_quad,"model_poisson"=fit_poisson,"model_quasi"=fit_quasi,"model_nb_lin"=fit_nb_lin,
                    "model_zinb" = fit_zinbh,"model_zinp" = fit_zinp,
                    "hnp_zinb" = my_hnp_zinbh,
                    "hnp_zinp" = my_hnp_zinp,
                    "hnp_nb_quad" = my_hnp_nb_quad,"hnp_nb_lin"=my_hnp_nb_lin,
                    "hnp_poisson" = my_hnp_poisson,
                    "hnp_quasi" =my_hnp_quasi,
                    
                    
                    "dist_nb_lin"=distance_nb_lin,
                    "dist_nb_quad" = distance_nb_quad,
                    "dist_poisson" = distance_poisson,
                    "dist_quasi" = distance_quasi,
                    "dist_zinbh" = distance_zinbh,"dist_zinp" = distance_zinp,
                    
                    
                    "dist_w_nb_lin"=distance_w_nb_lin,
                    "dist_w_nb_quad" = distance_w_nb_quad,
                    "dist_w_poisson" = distance_w_poisson,
                    "dist_w_quasi" = distance_w_quasi,
                    "dist_w_zinbh" = distance_w_zinbh,"dist_w_zinp" = distance_w_zinp,
                    
                    
                    
                    
                    "b_nb_quad"=b_nb_quad,"w_nb_quad"=w_nb_quad,"b_poisson"=b_poisson,"b_zinbh"=b_zinbh,"b_zinp"=b_zinp,
                    "w_poisson"=w_poisson,"b_quasi"=b_quasi,"w_quasi"=w_quasi,"w_zinbh"=w_zinbh,"w_zinp"=w_zinp,
                    "b_nb_lin"=b_nb_lin,"w_nb_lin"=w_nb_lin,
                    "aic_neg_lin"=aic_neg_lin,"aic_neg_quad"=aic_neg_quad,"aic_poisson"=aic_poisson,"aic_zinbh"=aic_zinbh,"aic_zinp"=aic_zinp,
                    "bic_neg_lin"=bic_neg_lin,"bic_neg_quad"=bic_neg_quad,"bic_poisson"=bic_poisson,"bic_zinbh"=bic_zinbh, "bic_zinp"=bic_zinp, 
                    "base_value_poisson" = base_value_poisson, "base_value_negbinlin" = base_value_negbinlin, "base_value_negbinquad" = base_value_negbinquad, " base_value_quasi" = base_value_quasi," base_value_zinbh" = base_value_zinbh,"base_value_zinp" = base_value_zinp,
                    "w_value_poisson" = w_value_poisson,  "w_value_negbinlin" = w_value_negbinlin, "w_value_negbinquad" = w_value_negbinquad, "w_value_quasi"= w_value_quasi , "w_value_zinbh"= w_value_zinbh ,"w_value_zinp"= w_value_zinp ,
                    "b_value_poisson" = b_value_poisson,  "b_value_negbinlin "= b_value_negbinlin,"b_value_negbinquad" = b_value_negbinquad,  "b_value_quasi" = b_value_quasi,"b_value_zinbh" = b_value_zinbh,"b_value_zinp" = b_value_zinp,
                    "r_poisson"= r_poisson, "r_negbinquad"= r_negbinquad, "r_negbinlin" = r_negbinlin, "r_quasi"= r_quasi,"r_zinbh"= r_zinbh,"r_zinp"= r_zinp,
                    "m_poisson"= m_poisson, "m_negbinquad"= m_negbinquad, "m_negbinlin" = m_negbinlin, "m_quasi"= m_quasi,"m_zinbh"= m_zinbh,"m_zinp"= m_zinp,
                    "nu_hat" =nu_hat
  )
  
  
  return(invisible(my_return))
}

# Set seed for reproducibility
set.seed(123)

# Run the simulation
my_simulation1 <- func_sim(n = 20, z = 1000)

# Save the results
save(my_simulation1, file = "sim_zinboh1_new.Rdata")
