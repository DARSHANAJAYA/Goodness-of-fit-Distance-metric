The Goodness of fit diagnoistic can divided into following parts:
1. Simulation files
   The simulation files are named as Poisson.R, Negbinh.R, Negbinlin.R, Zeroinflatedpoisson.R, Zeroinflatednegativebinomial.R, Underdispersed.R for simulating from poisson, Negative binomial with quadratic
   variance, Negative binomial with linear variance, zeroinflated poisson, zero inflated Negative binomial and quasi poisson respectively.

   For Negative binomial with quadratic variance for the conditions of mild and strong overdispersion , theta values are given as 2 and 0.142 respectively.
   For Negative binomial with linear variance for the conditions of mild and strong overdispersion, sigma values are given as 0.5 and 7 respectively.
   For zeroinflation values 0.2 and 0.6 are used for mild and strong zeroinflation and for mild and strong strong overdispersion 0.5 and 7 are used.
2. Data wrangling
   
datawrangling.R file does the datawrangling to formulate dataframes from the simulations Rdatafiles. For running this files you need to run the simulation accordingly and for L2 norm use j = 1:3 and for L1 norm use j = 46:48. Change file locations accordingly if simulated files are used.

3. Plots
   The plots.R code will give the final plot results. 
