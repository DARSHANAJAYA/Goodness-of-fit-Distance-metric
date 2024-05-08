rm(list= ls(all = TRUE))
set.seed(123)
## loading packages
library(hnp)
library(glmmTMB)
library(tidyverse)

## reading datasets
walleye <- read.csv("walleyePP.csv", stringsAsFactors = TRUE, header = TRUE) %>%
  mutate(YEAR = as.factor(YEAR))
walleye <- walleye %>%
  mutate(obs = factor(1:nrow(walleye)))



walleye12 <- walleye %>%
  filter(YEAR == "2012")
walleye12 <- walleye12 %>%
  mutate(obs = factor(1:nrow(walleye12)))

mean(walleye12$N)
var(walleye12$N)
summary(aov(N~AGE, data = walleye12))
walleye17 <- walleye %>%
  filter(YEAR == "2017")
walleye17 <- walleye17 %>%
  mutate(obs = factor(1:nrow(walleye17)))

fit1_b <- glm(N~AGE, family = poisson, data = walleye12)
anova(fit1_b, test = "Chisq")

####################################################################
## Model fitting
####################################################################

## Walleye full models
walleye_pois <- glm(N ~ AGE * YEAR, family = poisson, data = walleye)
walleye_qpois <- glm(N ~ AGE * YEAR, family = quasipoisson, data = walleye)
walleye_gp <- glmmTMB(N ~ AGE * YEAR, family = genpois, data = walleye)
walleye_nb1 <- glmmTMB(N ~ AGE * YEAR, family = nbinom1, data = walleye)
walleye_nb2 <- glm.nb(N ~ AGE * YEAR, data = walleye)
walleye_comp <- glmmTMB(N ~ AGE * YEAR, family = compois, data = walleye)
walleye_pn <- glmmTMB(N ~ AGE * YEAR + (1 | obs), family = poisson, data = walleye)

## Walleye 2012 models
walleye12_pois <- glm(N ~ AGE, family = poisson, data = walleye12)
walleye12_qpois <- glm(N ~ AGE, family = quasipoisson, data = walleye12)
walleye12_gp <- glmmTMB(N ~ AGE, family = genpois, data = walleye12)
walleye12_nb1 <- glmmTMB(N ~ AGE, family = nbinom1, data = walleye12)
walleye12_nb2 <- glm.nb(N ~ AGE, data = walleye12)
walleye12_comp <- glmmTMB(N ~ AGE, family = compois, data = walleye12)
walleye12_pn <- glmmTMB(N ~ AGE + (1 | obs), family = poisson, data = walleye12)

## Walleye 2017 models
walleye17_pois <- glm(N ~ AGE, family = poisson, data = walleye17)
walleye17_qpois <- glm(N ~ AGE, family = quasipoisson, data = walleye17)
walleye17_gp <- glmmTMB(N ~ AGE, family = genpois, data = walleye17)
walleye17_nb1 <- glmmTMB(N ~ AGE, family = nbinom1, data = walleye17)
walleye17_nb2 <- glm.nb(N ~ AGE, data = walleye17)
walleye17_comp <- glmmTMB(N ~ AGE, family = compois, data = walleye17)
walleye17_pn <- glmmTMB(N ~ AGE + (1 | obs), family = poisson, data = walleye17)

#########################################################################
## Computing and plotting the half-normal plots with a simulated envelope
#########################################################################

## Walleye full models - half-normal plots


dfun <- function(obj) {
  residuals(obj, type = "pearson")
}
sfun <- function(n, obj) {
  simulate(obj)[[1]]
}
ffun_gp <- function(response) {
  fit <- try(glmmTMB(response ~ AGE * YEAR, family = genpois, data = walleye), silent = TRUE)
  while(class(fit) == "try-error") {
    response2 <- sfun(1, walleye_gp)
    fit <- try(glmmTMB(response2 ~ AGE * YEAR, family = genpois, data = walleye), silent = TRUE)
  }
  return(fit)
}
ffun_nb1 <- function(response) {
  fit <- try(glmmTMB(response ~ AGE * YEAR, family = nbinom1, data = walleye), silent = TRUE)
  while(class(fit) == "try-error") {
    response2 <- sfun(1, walleye_gp)
    fit <- try(glmmTMB(response2 ~ AGE * YEAR, family = nbinom1, data = walleye), silent = TRUE)
  }
  return(fit)
}
ffun_comp <- function(response) {
  fit <- try(glmmTMB(response ~ AGE * YEAR, family = compois, data = walleye), silent = TRUE)
  while(class(fit) == "try-error") {
    response2 <- sfun(1, walleye_comp)
    fit <- try(glmmTMB(response2 ~ AGE * YEAR, family = compois, data = walleye), silent = TRUE)
  }
  return(fit)
}
ffun_pn <- function(response) {
  fit <- try(glmmTMB(response ~ AGE * YEAR + (1 | obs), family = poisson, data = walleye), silent = TRUE)
  while(class(fit) == "try-error") {
    response2 <- sfun(1, walleye_pn)
    fit <- try(glmmTMB(response2 ~ AGE * YEAR + (1 | obs), family = poisson, data = walleye), silent = TRUE)
  }
  return(fit)
}

set.seed(2020)
hnp_walleye_pois <- hnp(walleye_pois, plot = FALSE, resid.type = "pearson")
hnp_walleye_qpois <- hnp(walleye_qpois, plot = FALSE, resid.type = "pearson")
hnp_walleye_gp <- hnp(walleye_gp, plot = FALSE, newclass = TRUE,
                      diagfun = dfun, simfun = sfun, fitfun = ffun_gp)
hnp_walleye_nb1 <- hnp(walleye_nb1, plot = FALSE, newclass = TRUE,
                       diagfun = dfun, simfun = sfun, fitfun = ffun_nb1)
hnp_walleye_nb2 <- hnp(walleye_nb2, plot = FALSE, resid.type = "pearson")
hnp_walleye_comp <- hnp(walleye_comp, plot = FALSE, newclass = TRUE,
                        diagfun = dfun, simfun = sfun, fitfun = ffun_comp)
hnp_walleye_pn <- hnp(walleye_pn, plot = FALSE, newclass = TRUE,
                      diagfun = dfun, simfun = sfun, fitfun = ffun_pn)

hnp_walleye <- data_frame(residuals = c(hnp_walleye_pois$residuals,
                                        hnp_walleye_qpois$residuals,
                                        hnp_walleye_gp$residuals,
                                        hnp_walleye_nb1$residuals,
                                        hnp_walleye_nb2$residuals,
                                        hnp_walleye_comp$residuals,
                                        hnp_walleye_pn$residuals),
                          lower = c(hnp_walleye_pois$lower,
                                    hnp_walleye_qpois$lower,
                                    hnp_walleye_gp$lower,
                                    hnp_walleye_nb1$lower,
                                    hnp_walleye_nb2$lower,
                                    hnp_walleye_comp$lower,
                                    hnp_walleye_pn$lower),
                          median = c(hnp_walleye_pois$median,
                                     hnp_walleye_qpois$median,
                                     hnp_walleye_gp$median,
                                     hnp_walleye_nb1$median,
                                     hnp_walleye_nb2$median,
                                     hnp_walleye_comp$median,
                                     hnp_walleye_pn$median),
                          upper = c(hnp_walleye_pois$upper,
                                    hnp_walleye_qpois$upper,
                                    hnp_walleye_gp$upper,
                                    hnp_walleye_nb1$upper,
                                    hnp_walleye_nb2$upper,
                                    hnp_walleye_comp$upper,
                                    hnp_walleye_pn$upper),
                          x = c(hnp_walleye_pois$x,
                                hnp_walleye_qpois$x,
                                hnp_walleye_gp$x,
                                hnp_walleye_nb1$x,
                                hnp_walleye_nb2$x,
                                hnp_walleye_comp$x,
                                hnp_walleye_pn$x),
                          model = factor(rep(c("Poisson",
                                               "Quasi-Poisson",
                                               "Generalized Poisson",
                                               "Negative binomial type 1",
                                               "Negative binomial type 2",
                                               "Mean-parameterized Conway-Maxwell-Poisson",
                                               "Poisson-normal"), each = 99),
                                         levels = c("Poisson",
                                                    "Quasi-Poisson",
                                                    "Negative binomial type 1",
                                                    "Negative binomial type 2",
                                                    "Mean-parameterized Conway-Maxwell-Poisson",
                                                    "Generalized Poisson",
                                                    "Poisson-normal")))

plot_walleye <- hnp_walleye %>%
  ggplot(aes(x = x, y = residuals)) +
  theme_bw() +
  facet_wrap(~ model, scales = "free") +
  geom_point(cex = 1, pch = 16, alpha = .75) +
  geom_line(aes(y = median),
            lty = 2, lwd = .2) +
  geom_line(aes(y = upper),
            lty = 1, lwd = .2) +
  geom_line(aes(y = lower),
            lty = 1, lwd = .2) +
  geom_ribbon(aes(ymin = lower, ymax = upper),
              fill = "gray", alpha = .2) +
  xlab("Half-normal scores") +
  ylab("Pearson residuals") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  ggtitle("Full model")

## Walleye 2012 models - half-normal plots
dfun <- function(obj) {
  residuals(obj, type = "pearson")
}
sfun <- function(n, obj) {
  simulate(obj)[[1]]
}
ffun_nb1 <- function(response) {
  fit <- try(glmmTMB(response ~ AGE, family = nbinom1, data = walleye12), silent = TRUE)
  while(class(fit) == "try-error") {
    response2 <- sfun(1, walleye12_gp)
    fit <- try(glmmTMB(response2 ~ AGE, family = nbinom1, data = walleye12), silent = TRUE)
  }
  return(fit)
}
ffun_gp <- function(response) {
  fit <- try(glmmTMB(response ~ AGE, family = genpois, data = walleye12), silent = TRUE)
  while(class(fit) == "try-error") {
    response2 <- sfun(1, walleye12_gp)
    fit <- try(glmmTMB(response2 ~ AGE, family = genpois, data = walleye12), silent = TRUE)
  }
  return(fit)
}
ffun_comp <- function(response) {
  fit <- try(glmmTMB(response ~ AGE, family = compois, data = walleye12), silent = TRUE)
  while(class(fit) == "try-error") {
    response2 <- sfun(1, walleye12_comp)
    fit <- try(glmmTMB(response2 ~ AGE, family = compois, data = walleye12), silent = TRUE)
  }
  return(fit)
}
ffun_pn <- function(response) {
  fit <- try(glmmTMB(response ~ AGE + (1 | obs), family = poisson, data = walleye12), silent = TRUE)
  while(class(fit) == "try-error") {
    response2 <- sfun(1, walleye12_pn)
    fit <- try(glmmTMB(response2 ~ AGE + (1 | obs), family = poisson, data = walleye12), silent = TRUE)
  }
  return(fit)
}

set.seed(2020)
hnp_walleye12_pois <- hnp(walleye12_pois, plot = FALSE, resid.type = "pearson")
hnp_walleye12_qpois <- hnp(walleye12_qpois, plot = FALSE, resid.type = "pearson")
hnp_walleye12_gp <- hnp(walleye12_gp, plot = FALSE, newclass = TRUE,
                        diagfun = dfun, simfun = sfun, fitfun = ffun_gp)
hnp_walleye12_nb1 <- hnp(walleye12_nb1, plot = FALSE, newclass = TRUE,
                         diagfun = dfun, simfun = sfun, fitfun = ffun_nb1)



hnp_walleye12_nb2 <- hnp(walleye12_nb2, plot = FALSE, resid.type = "pearson")
hnp_walleye12_comp <- hnp(walleye12_comp, plot = FALSE, newclass = TRUE,
                          diagfun = dfun, simfun = sfun, fitfun = ffun_comp)
hnp_walleye12_pn <- hnp(walleye12_pn, plot = FALSE, newclass = TRUE,
                        diagfun = dfun, simfun = sfun, fitfun = ffun_pn)

hnp_walleye12 <- data_frame(residuals = c(hnp_walleye12_pois$residuals,
                                          hnp_walleye12_qpois$residuals,
                                          hnp_walleye12_gp$residuals,
                                          hnp_walleye12_nb1$residuals,
                                          hnp_walleye12_nb2$residuals,
                                          hnp_walleye12_comp$residuals,
                                          hnp_walleye12_pn$residuals),
                            lower = c(hnp_walleye12_pois$lower,
                                      hnp_walleye12_qpois$lower,
                                      hnp_walleye12_gp$lower,
                                      hnp_walleye12_nb1$lower,
                                      hnp_walleye12_nb2$lower,
                                      hnp_walleye12_comp$lower,
                                      hnp_walleye12_pn$lower),
                            median = c(hnp_walleye12_pois$median,
                                       hnp_walleye12_qpois$median,
                                       hnp_walleye12_gp$median,
                                       hnp_walleye12_nb1$median,
                                       hnp_walleye12_nb2$median,
                                       hnp_walleye12_comp$median,
                                       hnp_walleye12_pn$median),
                            upper = c(hnp_walleye12_pois$upper,
                                      hnp_walleye12_qpois$upper,
                                      hnp_walleye12_gp$upper,
                                      hnp_walleye12_nb1$upper,
                                      hnp_walleye12_nb2$upper,
                                      hnp_walleye12_comp$upper,
                                      hnp_walleye12_pn$upper),
                            x = c(hnp_walleye12_pois$x,
                                  hnp_walleye12_qpois$x,
                                  hnp_walleye12_gp$x,
                                  hnp_walleye12_nb1$x,
                                  hnp_walleye12_nb2$x,
                                  hnp_walleye12_comp$x,
                                  hnp_walleye12_pn$x),
                            model = factor(rep(c("Poisson",
                                                 "Quasi-Poisson",
                                                 "Generalized Poisson",
                                                 "Negative binomial type 1",
                                                 "Negative binomial type 2",
                                                 "Mean-parameterized Conway-Maxwell-Poisson",
                                                 "Poisson-normal"), each = 44),
                                           levels = c("Poisson",
                                                      "Quasi-Poisson",
                                                      "Negative binomial type 1",
                                                      "Negative binomial type 2",
                                                      "Mean-parameterized Conway-Maxwell-Poisson",
                                                      "Generalized Poisson",
                                                      "Poisson-normal")))

plot_walleye12 <- hnp_walleye12 %>%
  ggplot(aes(x = x, y = residuals)) +
  theme_bw() +
  facet_wrap(~ model, scales = "free") +
  geom_point(cex = 1, pch = 16, alpha = .75) +
  geom_line(aes(y = median),
            lty = 2, lwd = .2) +
  geom_line(aes(y = upper),
            lty = 1, lwd = .2) +
  geom_line(aes(y = lower),
            lty = 1, lwd = .2) +
  geom_ribbon(aes(ymin = lower, ymax = upper),
              fill = "gray", alpha = .2) +
  xlab("Half-normal scores") +
  ylab("Pearson residuals") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  ggtitle("Walleye 2012")

## Walleye 2017 models - half-normal plots
dfun <- function(obj) {
  residuals(obj, type = "pearson")
}
sfun <- function(n, obj) {
  simulate(obj)[[1]]
}
ffun_nb1 <- function(response) {
  fit <- try(glmmTMB(response ~ AGE, family = nbinom1, data = walleye17), silent = TRUE)
  while(class(fit) == "try-error") {
    response2 <- sfun(1, walleye17_gp)
    fit <- try(glmmTMB(response2 ~ AGE, family = nbinom1, data = walleye17), silent = TRUE)
  }
  return(fit)
}
ffun_gp <- function(response) {
  fit <- try(glmmTMB(response ~ AGE, family = genpois, data = walleye17), silent = TRUE)
  while(class(fit) == "try-error") {
    response2 <- sfun(1, walleye17_gp)
    fit <- try(glmmTMB(response2 ~ AGE, family = genpois, data = walleye17), silent = TRUE)
  }
  return(fit)
}
ffun_comp <- function(response) {
  fit <- try(glmmTMB(response ~ AGE, family = compois, data = walleye17), silent = TRUE)
  while(class(fit) == "try-error") {
    response2 <- sfun(1, walleye17_comp)
    fit <- try(glmmTMB(response2 ~ AGE, family = compois, data = walleye17), silent = TRUE)
  }
  return(fit)
}
ffun_pn <- function(response) {
  fit <- try(glmmTMB(response ~ AGE + (1 | obs), family = poisson, data = walleye17), silent = TRUE)
  while(class(fit) == "try-error") {
    response2 <- sfun(1, walleye17_pn)
    fit <- try(glmmTMB(response2 ~ AGE + (1 | obs), family = poisson, data = walleye17), silent = TRUE)
  }
  return(fit)
}

set.seed(2020)
hnp_walleye17_pois <- hnp(walleye17_pois, plot = FALSE, resid.type = "pearson")
hnp_walleye17_qpois <- hnp(walleye17_qpois, plot = FALSE, resid.type = "pearson")
hnp_walleye17_gp <- hnp(walleye17_gp, plot = FALSE, newclass = TRUE,
                        diagfun = dfun, simfun = sfun, fitfun = ffun_gp)
hnp_walleye17_nb1 <- hnp(walleye17_nb1, plot = FALSE, newclass = TRUE,
                         diagfun = dfun, simfun = sfun, fitfun = ffun_nb1)
 hnp_walleye17_nb2 <- hnp(walleye17_nb2, plot = FALSE, resid.type = "pearson")
hnp_walleye17_comp <- hnp(walleye17_comp, plot = FALSE, newclass = TRUE,
                          diagfun = dfun, simfun = sfun, fitfun = ffun_comp)
hnp_walleye17_pn <- hnp(walleye17_pn, plot = FALSE, newclass = TRUE,
                        diagfun = dfun, simfun = sfun, fitfun = ffun_pn)

hnp_walleye17 <- data_frame(residuals = c(hnp_walleye17_pois$residuals,
                                          hnp_walleye17_qpois$residuals,
                                          hnp_walleye17_gp$residuals,
                                          hnp_walleye17_nb1$residuals,
                                          hnp_walleye17_nb2$residuals,
                                          hnp_walleye17_comp$residuals,
                                          hnp_walleye17_pn$residuals),
                            lower = c(hnp_walleye17_pois$lower,
                                      hnp_walleye17_qpois$lower,
                                      hnp_walleye17_gp$lower,
                                      hnp_walleye17_nb1$lower,
                                      hnp_walleye17_nb2$lower,
                                      hnp_walleye17_comp$lower,
                                      hnp_walleye17_pn$lower),
                            median = c(hnp_walleye17_pois$median,
                                       hnp_walleye17_qpois$median,
                                       hnp_walleye17_gp$median,
                                       hnp_walleye17_nb1$median,
                                       hnp_walleye17_nb2$median,
                                       hnp_walleye17_comp$median,
                                       hnp_walleye17_pn$median),
                            upper = c(hnp_walleye17_pois$upper,
                                      hnp_walleye17_qpois$upper,
                                      hnp_walleye17_gp$upper,
                                      hnp_walleye17_nb1$upper,
                                      hnp_walleye17_nb2$upper,
                                      hnp_walleye17_comp$upper,
                                      hnp_walleye17_pn$upper),
                            x = c(hnp_walleye17_pois$x,
                                  hnp_walleye17_qpois$x,
                                  hnp_walleye17_gp$x,
                                  hnp_walleye17_nb1$x,
                                  hnp_walleye17_nb2$x,
                                  hnp_walleye17_comp$x,
                                  hnp_walleye17_pn$x),
                            model = factor(rep(c("Poisson",
                                                 "Quasi-Poisson",
                                                 "Generalized Poisson",
                                                 "Negative binomial type 1",
                                                 "Negative binomial type 2",
                                                 "Mean-parameterized Conway-Maxwell-Poisson",
                                                 "Poisson-normal"), each = 55),
                                           levels = c("Poisson",
                                                      "Quasi-Poisson",
                                                      "Negative binomial type 1",
                                                      "Negative binomial type 2",
                                                      "Mean-parameterized Conway-Maxwell-Poisson",
                                                      "Generalized Poisson",
                                                      "Poisson-normal")))

plot_walleye17 <- hnp_walleye17 %>%
  ggplot(aes(x = x, y = residuals)) +
  theme_bw() +
  facet_wrap(~ model, scales = "free") +
  geom_point(cex = 1, pch = 16, alpha = .75) +
  geom_line(aes(y = median),
            lty = 2, lwd = .2) +
  geom_line(aes(y = upper),
            lty = 1, lwd = .2) +
  geom_line(aes(y = lower),
            lty = 1, lwd = .2) +
  geom_ribbon(aes(ymin = lower, ymax = upper),
              fill = "gray", alpha = .2) +
  xlab("Half-normal scores") +
  ylab("Pearson residuals") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  ggtitle("Walleye-2017")

######## Distance metric function

my_distance<-function(hnp_obj) {
  r <-hnp_obj$residuals
  m<-hnp_obj$median
  
  
  # --- old code
  #distance<-sum(((r-m)^2/(f))*(g)^Indicator)
  
  # --- new code
  distance<-sum((r-m)^2)
  distance_1 <- sum(abs(r-m))
  
  my_list<-list(distance,distance_1)
  
  return(my_list)
}


########Distance metric for walleye full models 

d_pois <- my_distance(hnp_walleye_pois)
d_qpois <- my_distance(hnp_walleye_qpois)
d_gp <- my_distance(hnp_walleye_gp)
# d_nb1 <- my_distance(hnp_walleye_nb1)
d_nb2 <- my_distance(hnp_walleye_nb2)
d_comp <- my_distance(hnp_walleye_comp)
d_pn <- my_distance(hnp_walleye_pn)


#######Distance metric for walleye 2012 models

d_pois_12 <- my_distance(hnp_walleye12_pois)
d_qpois_12 <- my_distance(hnp_walleye12_qpois)
d_gp_12 <- my_distance(hnp_walleye12_gp)
d_nb1_12  <- my_distance(hnp_walleye12_nb1)
d_nb2_12  <- my_distance(hnp_walleye12_nb2)
d_comp_12  <- my_distance(hnp_walleye12_comp)
d_pn_12  <- my_distance(hnp_walleye12_pn)


######Distance metric for walleye 2017 models 

d_pois_17 <- my_distance(hnp_walleye17_pois)
d_qpois_17 <- my_distance(hnp_walleye17_qpois)
d_gp_17 <- my_distance(hnp_walleye17_gp)
d_nb1_17  <- my_distance(hnp_walleye17_nb1)
d_nb2_17  <- my_distance(hnp_walleye17_nb2)
d_comp_17  <- my_distance(hnp_walleye17_comp)
d_pn_17  <- my_distance(hnp_walleye17_pn)




####################New code for 100 simulations of the above code 


rm(list= ls(all = TRUE))
set.seed(123)
## loading packages
library(hnp)
library(glmmTMB)
library(tidyverse)

## reading datasets
walleye <- read.csv("walleyePP.csv", stringsAsFactors = TRUE, header = TRUE) %>%
  mutate(YEAR = as.factor(YEAR))
walleye <- walleye %>%
  mutate(obs = factor(1:nrow(walleye)))



walleye12 <- walleye %>%
  filter(YEAR == "2012")
walleye12 <- walleye12 %>%
  mutate(obs = factor(1:nrow(walleye12)))
walleye17 <- walleye %>%
  filter(YEAR == "2017")
walleye17 <- walleye17 %>%
  mutate(obs = factor(1:nrow(walleye17)))

ggplot(walleye12, aes(y = N, x  = AGE)) + 
  geom_point()
walleye12 %>%
  ggplot(aes(x = AGE, y = N)) +
  theme_bw() +
  geom_point() +
  geom_line() 

####################################################################
## Model fitting
####################################################################

## Walleye full models
walleye_pois <- glm(N ~ AGE * YEAR, family = poisson, data = walleye)
walleye_qpois <- glm(N ~ AGE * YEAR, family = quasipoisson, data = walleye)
walleye_nb2 <- glm.nb(N ~ AGE * YEAR, data = walleye)
walleye_zip <- zeroinfl(N ~ AGE, dist = "poisson", data = walleye)
walleye_zinb <- zeroinfl(N ~ AGE, dist = "negbin", data = walleye)

## Walleye 2012 models
walleye12_pois <- glm(N ~ AGE, family = poisson, data = walleye12)
walleye12_qpois <- glm(N ~ AGE, family = quasipoisson, data = walleye12)
walleye12_nblin <- gamlss(N ~ AGE, family = NBII, data = walleye12)
walleye12_nb2 <- glm.nb(N ~ AGE, data = walleye12)
walleye12_zip <- zeroinfl(N ~ AGE , dist = "poisson", data = walleye12)
walleye12_zinb <- zeroinfl(N ~ AGE, dist = "negbin", data = walleye12)
BIC(walleye12_pois)
BIC(walleye12_qpois)
BIC(walleye12_nblin)
BIC(walleye12_nb2)
BIC(walleye12_zip)
BIC(walleye12_zinb)

par(mfrow = c(2, 3))

hnp12_p <- hnp(walleye12_pois, main = "Poisson",xlab = "Half normal scores", resid.type = "pearson")
hnp12_qp <- hnp(walleye12_qpois, main = "Quasi", xlab = "Half normal scores", resid.type = "pearson")


####NB LIN


dfun <- function(obj) {
  r <- obj$y - obj$mu.fv
  v <- obj$mu.fv * (1 + obj$sigma.fv)
  rp <- r / sqrt(v)
  return(rp)
}
sfun <- function(n, obj) {
  y <- rNBII(length(obj$mu.fv), mu = obj$mu.fv, sigma = obj$sigma.fv)
  return(y)
}
ffun <- function(resp) {
  walleye12$resp <- resp
  fit <- gamlss(resp ~ AGE, family = NBII, data = walleye12)
  return(fit)
}

hnp12_nblin <- hnp(walleye12_nblin,newclass = TRUE,diagfun = dfun, fitfun = ffun, 
                   simfun = sfun, main = "NB-lin",xlab = "Half normal scores",  resid.type = "pearson")

hnp12_nb <- hnp(walleye12_nb2, main = "NB-Quad",xlab = "Half normal scores",  resid.type = "pearson")
hnp12_zip <- hnp(walleye12_zip, main = "ZIP", xlab = "Half normal scores", resid.type = "pearson")
hnp12_zinb <- hnp(walleye12_zinb, main = "ZINB",xlab = "Half normal scores", resid.type = "pearson")





## Walleye 2017 models
walleye17_pois <- glm(N ~ AGE, family = poisson, data = walleye17)
walleye17_qpois <- glm(N ~ AGE, family = quasipoisson, data = walleye17)
walleye17_nb2 <- glm.nb(N ~ AGE, data = walleye17)
walleye17_zip <- zeroinfl(N ~ AGE, dist = "poisson", data = walleye17)
walleye17_zinb <- zeroinfl(N ~ AGE, dist = "negbin", data = walleye17)

my_distance<-function(hnp_obj) {
  r <-hnp_obj$residuals
  m<-hnp_obj$median
  
  
  # --- old code
  #distance<-sum(((r-m)^2/(f))*(g)^Indicator)
  
  # --- new code
  distance<-sum((r-m)^2)
  distance_1 <- sum(abs(r-m))
  
  my_list<-list(distance,distance_1)
  
  return(my_list)
  
}

distance_p_l2 <- distance_nb_l2 <- distance_qp_l2 <- distance_zinb_l2 <- distance_zip_l2 <- list()

distance_p_l1 <- distance_nb_l1 <- distance_qp_l1 <- distance_zinb_l1 <- distance_zip_l1 <- list()

distance12_p_l2 <- distance12_nb_l2 <- distance12_qp_l2 <- distance12_zinb_l2 <- distance12_zip_l2 <- list()

distance12_p_l1 <- distance12_nb_l1 <- distance12_qp_l1 <- distance12_zinb_l1 <- distance12_zip_l1 <- list()

distance17_p_l2 <- distance17_nb_l2 <- distance17_qp_l2 <- distance17_zinb_l2 <- distance17_zip_l2 <- list()

distance17_p_l1 <- distance17_nb_l1 <- distance17_qp_l1 <- distance17_zinb_l1 <- distance17_zip_l1 <- list()


for(i in 1:100){
  # hnp_p <- hnp(walleye_pois, main = "Poisson",xlab = "Half normal scores", plot = FALSE, resid.type = "pearson")
  # hnp_qp <- hnp(walleye_qpois, main = "Quasi", xlab = "Half normal scores", plot = FALSE, resid.type = "pearson")
  # hnp_zip <- hnp(walleye_zip, main = "ZIP", xlab = "Half normal scores", plot = FALSE, resid.type = "pearson")
  # hnp_zinb <- hnp(walleye_zinb, main = "ZINB",xlab = "Half normal scores", plot = FALSE, resid.type = "pearson")
  # hnp_nb <- hnp(walleye_nb2, main = "NB-Quad",xlab = "Half normal scores", plot = FALSE, resid.type = "pearson")
  # 
  hnp12_p <- hnp(walleye12_pois, main = "Poisson",xlab = "Half normal scores", plot = FALSE, resid.type = "pearson")
  hnp12_qp <- hnp(walleye12_qpois, main = "Quasi", xlab = "Half normal scores", plot = FALSE, resid.type = "pearson")
  
  dfun <- function(obj) {
    r <- obj$y - obj$mu.fv
    v <- obj$mu.fv * (1 + obj$sigma.fv)
    rp <- r / sqrt(v)
    return(rp)
  }
  sfun <- function(n, obj) {
    y <- rNBII(length(obj$mu.fv), mu = obj$mu.fv, sigma = obj$sigma.fv)
    return(y)
  }
  ffun <- function(resp) {
    walleye12$resp <- resp
    fit <- gamlss(resp ~ AGE, family = NBII, data = walleye12)
    return(fit)
  }
  
  # hnp12_nblin <- hnp(fit_, newclass = TRUE,diagfun = dfun, fitfun = ffun, simfun = sfun)
  hnp12_zip <- hnp(walleye12_zip, main = "ZIP", xlab = "Half normal scores", plot = FALSE, resid.type = "pearson")
  hnp12_zinb <- hnp(walleye12_zinb, main = "ZINB",xlab = "Half normal scores", plot = FALSE, resid.type = "pearson")
  hnp12_nb <- hnp(walleye12_nb2, main = "NB-Quad",xlab = "Half normal scores", plot = FALSE, resid.type = "pearson")
  
  
  # hnp17_p <- hnp(walleye17_pois, main = "Poisson",xlab = "Half normal scores", plot = FALSE, resid.type = "pearson")
  # hnp17_qp <- hnp(walleye17_qpois, main = "Quasi", xlab = "Half normal scores", plot = FALSE, resid.type = "pearson")
  # hnp17_zip <- hnp(walleye17_zip, main = "ZIP", xlab = "Half normal scores", plot = FALSE, resid.type = "pearson")
  # hnp17_zinb <- hnp(walleye17_zinb, main = "ZINB",xlab = "Half normal scores", plot = FALSE, resid.type = "pearson")
  # hnp17_nb <- hnp(walleye17_nb2, main = "NB-Quad",xlab = "Half normal scores", plot = FALSE, resid.type = "pearson")
  
  
  # distance_p_l2[i] <- my_distance(hnp_p)[1]
  # distance_nb_l2[i] <- my_distance(hnp_nb)[1]
  # distance_qp_l2[i] <- my_distance(hnp_qp)[1]
  # distance_zinb_l2[i] <- my_distance(hnp_zinb)[1]
  # distance_zip_l2[i] <- my_distance(hnp_zip)[1]
  # 
  # distance_p_l1[i] <- my_distance(hnp_p)[2]
  # distance_nb_l1[i] <- my_distance(hnp_nb)[2]
  # distance_qp_l1[i] <- my_distance(hnp_qp)[2]
  # distance_zinb_l1[i] <- my_distance(hnp_zinb)[2]
  # distance_zip_l1[i] <- my_distance(hnp_zip)[2]
  # 
  distance12_p_l2[i] <- my_distance(hnp12_p)[1]
  distance12_nb_l2[i] <- my_distance(hnp12_nb)[1]
  distance12_qp_l2[i] <- my_distance(hnp12_qp)[1]
  distance12_zinb_l2[i] <- my_distance(hnp12_zinb)[1]
  distance12_zip_l2[i] <- my_distance(hnp12_zip)[1]
  
  distance12_p_l1[i] <- my_distance(hnp12_p)[2]
  distance12_nb_l1[i] <- my_distance(hnp12_nb)[2]
  distance12_qp_l1[i] <- my_distance(hnp12_qp)[2]
  distance12_zinb_l1[i] <- my_distance(hnp12_zinb)[2]
  distance12_zip_l1[i] <- my_distance(hnp12_zip)[2]
  
  # distance17_p_l2[i] <- my_distance(hnp17_p)[1]
  # distance17_nb_l2[i] <- my_distance(hnp17_nb)[1]
  # distance17_qp_l2[i] <- my_distance(hnp17_qp)[1]
  # distance17_zinb_l2[i] <- my_distance(hnp17_zinb)[1]
  # distance17_zip_l2[i] <- my_distance(hnp17_zip)[1]
  # 
  # distance17_p_l1[i] <- my_distance(hnp17_p)[2]
  # distance17_nb_l1[i] <- my_distance(hnp17_nb)[2]
  # distance17_qp_l1[i] <- my_distance(hnp17_qp)[2]
  # distance17_zinb_l1[i] <- my_distance(hnp17_zinb)[2]
  # distance17_zip_l1[i] <- my_distance(hnp17_zip)[2]
  
}

distance_l2 <- tibble(distance12_p_l2, distance12_nb_l2, distance12_qp_l2, distance12_zinb_l2, distance12_zip_l2)
distance_l1 <- tibble(distance12_p_l1, distance12_nb_l1, distance12_qp_l1, distance12_zinb_l1, distance12_zip_l1)
# distance_l2 <- as.data.frame(distance_l2)
summary(unlist(distance_l2$distance12_p_l2))
IQR(unlist(distance_l2$distance12_p_l2))
sd(unlist(distance_l2$distance12_p_l2))
summary(unlist(distance_l2$distance12_nb_l2))
IQR(unlist(distance_l2$distance12_nb_l2))
sd(unlist(distance_l2$distance12_nb_l2))
summary(unlist(distance_l2$distance12_qp_l2))
IQR(unlist(distance_l2$distance12_qp_l2))
sd(unlist(distance_l2$distance12_qp_l2))
summary(unlist(distance_l2$distance12_zinb_l2))
IQR(unlist(distance_l2$distance12_zinb_l2))
sd(unlist(distance_l2$distance12_zinb_l2))
summary(unlist(distance_l2$distance12_zip_l2))
IQR(unlist(distance_l2$distance12_zip_l2))
sd(unlist(distance_l2$distance12_zip_l2))

###This is l1 values
summary(unlist(distance_l1$distance12_p_l1))

IQR(unlist(distance_l1$distance12_p_l1))
sd(unlist(distance_l1$distance12_p_l1))
summary(unlist(distance_l1$distance12_nb_l1))
IQR(unlist(distance_l1$distance12_nb_l1))
sd(unlist(distance_l1$distance12_nb_l1))
summary(unlist(distance_l1$distance12_qp_l1))
IQR(unlist(distance_l1$distance12_qp_l1))
sd(unlist(distance_l1$distance12_qp_l1))
summary(unlist(distance_l1$distance12_zinb_l1))
IQR(unlist(distance_l1$distance12_zinb_l1))
sd(unlist(distance_l1$distance12_zinb_l1))
summary(unlist(distance_l1$distance12_zip_l1))
IQR(unlist(distance_l1$distance12_zip_l1))
sd(unlist(distance_l1$distance12_zip_l1))

####TRYING from rafaels code 
  
  library(mvabund)
  library(gamlss)
  library(hnp)
  
  ## reading datasets
  walleye <- read.csv("walleyePP.csv", stringsAsFactors = TRUE, header = TRUE) %>%
    mutate(YEAR = as.factor(YEAR))
  walleye <- walleye %>%
    mutate(obs = factor(1:nrow(walleye)))

  
  walleye12 <- walleye %>%
    filter(YEAR == "2012")
  walleye12 <- walleye12 %>%
    mutate(obs = factor(1:nrow(walleye12)))
  
  fit <- glm(N ~ AGE, family = quasipoisson, data = walleye12)
  hnp(fit, resid.type = "pearson")
  
  fit2 <- gamlss(N ~ AGE, family = NBII, data = walleye12)
  summary(fit2)
  BIC(fit2)
  dfun <- function(obj) {
    r <- obj$y - obj$mu.fv
    v <- obj$mu.fv * (1 + obj$sigma.fv)
    rp <- r / sqrt(v)
    return(rp)
  }
  sfun <- function(n, obj) {
    y <- rNBII(length(obj$mu.fv), mu = obj$mu.fv, sigma = obj$sigma.fv)
    return(y)
  }
  ffun <- function(resp) {
    walleye12$resp <- resp
    fit <- gamlss(resp ~ AGE, family = NBII, data = walleye12)
    return(fit)
  }
  
  hnp(fit2, newclass = TRUE,
      diagfun = dfun, fitfun = ffun, simfun = sfun)
  
  
  
  ############RUNNING THE SIMULATION SEPERATELY FOR NB LIN 
  
  
  ########TRYING IT SEPERATELY WITH THE CODE 
  
  
  ## reading datasets
  # walleye <- read.csv("walleyePP.csv", stringsAsFactors = TRUE, header = TRUE) %>%
  #   mutate(YEAR = as.factor(YEAR))
  # walleye <- walleye %>%
  #   mutate(obs = factor(1:nrow(walleye)))
  
  
  walleye12 <- walleye %>%
    filter(YEAR == "2012")
  walleye12 <- walleye12 %>%
    mutate(obs = factor(1:nrow(walleye12)))
  
  # fit <- glm(y ~ x, family = quasipoisson, data = d)
  # hnp(fit, resid.type = "pearson")
  
  fit2 <- gamlss(N ~ AGE, family = NBII, data = walleye12)
  summary(fit2)
  distance_fit2_l2 <- distance_fit2_l1 <- list()
  
  my_distance<-function(hnp_obj) {
    r <-hnp_obj$residuals
    m<-hnp_obj$median
    
    
    # --- old code
    #distance<-sum(((r-m)^2/(f))*(g)^Indicator)
    
    # --- new code
    distance<-sum((r-m)^2)
    distance_1 <- sum(abs(r-m))
    
    my_list<-list(distance,distance_1)
    
    return(my_list)
    
  }
  
  for(i in 1:100){
    dfun <- function(obj) {
      r <- obj$y - obj$mu.fv
      v <- obj$mu.fv * (1 + obj$sigma.fv)
      rp <- r / sqrt(v)
      return(rp)
    }
    sfun <- function(n, obj) {
      y <- rNBII(length(obj$mu.fv), mu = obj$mu.fv, sigma = obj$sigma.fv)
      return(y)
    }
    ffun <- function(resp) {
      walleye12$resp <- resp
      fit <- gamlss(resp ~ AGE, family = NBII, data = walleye12)
      return(fit)
    }
    
    
    hnp_fit2 <- hnp(fit2, newclass = TRUE,
                    diagfun = dfun, fitfun = ffun, simfun = sfun, plot = FALSE)
    
    distance_fit2_l2[i] <- my_distance(hnp_fit2)[1]
    distance_fit2_l1[i] <- my_distance(hnp_fit2)[2]
    
  }
  
  ###########exponentiating the ZINB coeffecients 
  
  expCoef <- exp(coef((walleye12_zinb)))
  expCoef <- matrix(expCoef, ncol = 2)
  rownames(expCoef) <- names(coef(walleye12_zinb))
  colnames(expCoef) <- c("Count_model","Zero_inflation_model")
  expCoef

