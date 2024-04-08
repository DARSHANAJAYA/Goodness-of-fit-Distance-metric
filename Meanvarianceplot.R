# Load necessary library
require(gamlss)

# Set parameters
n <- 100 # Sample size
beta_coef <- c(3, 2) # Coefficients for the linear predictor
theta <- 5 # Dispersion parameter

# Generate data for negative binomial distribution
x1 <- rnorm(n) # Generate predictor variable
mu <- exp(beta_coef[1] + beta_coef[2]*x1) # Calculate the mean using the linear predictor
y1 <- rNBI(n = n, mu = mu, sigma = 1/theta) # Generate response variable following negative binomial distribution

# Plot generated data
plot(x1, y1)

# Group predictor variable into intervals
x1_cut <- cut(x1, breaks = seq(-3, 3, length = 20))

# Calculate means and variances for each group
means <- tapply(y1, x1_cut, mean)
vars <- tapply(y1, x1_cut, var)

# Plot means against variances
plot(means, vars)

# Load data from saved file
load("/Volumes/Darshana_PhD/project/updated/L2/l2_negbinh/sim_negbinh3_new.Rdata")
y1 <- my_simulation1$model_poisson[[1]]$y # Get response variable from loaded data
x1 <- model.matrix(my_simulation1$model_poisson[[1]])[,2] # Get predictor variable from loaded data

# Plot loaded data
plot(x1, y1)

# Group predictor variable into intervals
x1_cut <- cut(x1, breaks = seq(-3, 3, length = 20))

# Calculate means and variances for each group
means <- tapply(y1, x1_cut, mean)
vars <- tapply(y1, x1_cut, var)

# Plot means against variances
plot(means, vars)
