###clearing the workspace 
rm(list = ls(all.names = TRUE))
###setting the seed for reproducibility 
set.seed(123)
#loading the libraries 
library(tidyverse)
library(gridExtra)
library(dplyr)
library(ggplot2)
#loading the data 
spider_df <- read_csv("https://raw.githubusercontent.com/rafamoral/courses/main/intro_glm/data/spider.csv")

colnames(spider_df) <- c("counts", "soil dry mass")
##exploratory plot 
ggplot(spider_df, aes(y = counts, x = `soil dry mass`)) + 
  geom_point()

summary(spider_df)
###fitting the models and calculating the corresponding the hnp objects 
par(mfrow = c(2, 3))
fit_p <- glm(count ~ `soil dry mass`,
             family = poisson,
             data = spider_df)
hnp_p <- hnp(fit_p, main = "Poisson",xlab = "Half normal scores", resid.type = "pearson")

fit_qp <- glm(count ~ soil_dry_mass,
              family = quasipoisson,
              data = spider_df)
hnp_qp <- hnp(fit_qp, main = "Quasi", xlab = "Half normal scores", resid.type = "pearson")



fit_nblin <- gamlss(count ~ soil_dry_mass, family = NBII, data = spider_df)
### external functions to calculate the 
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
  spider_df$resp <- resp
  fit <- gamlss(resp ~ soil_dry_mass, family = NBII, data = spider_df)
  return(fit)
}

hnp(fit_nblin, newclass = TRUE,
    diagfun = dfun, fitfun = ffun, simfun = sfun, main = "NB-lin",xlab = "Half normal scores",  resid.type = "pearson")


fit_nb <- glm.nb(count ~ soil_dry_mass,
                 data = spider_df)
hnp_nb <- hnp(fit_nb, main = "NB-Quad",xlab = "Half normal scores",  resid.type = "pearson")

library(pscl)
fit_zip <- zeroinfl(count ~ soil_dry_mass,
                    dist = "poisson",
                    data = spider_df)
hnp_zip <- hnp(fit_zip, main = "ZIP", xlab = "Half normal scores",  resid.type = "pearson")


fit_zinb <- zeroinfl(count ~ soil_dry_mass,
                     dist = "negbin",
                     data = spider_df)
hnp_zinb <- hnp(fit_zinb, main = "ZINB",xlab = "Half normal scores", resid.type = "pearson")


# 
# 
# fit_nb_lin <-gamlss(count ~ soil_dry_mass, family = NBII,data = spider_df)#nblin
# d.fun <- function(obj) {
#   r <- obj$y - obj$mu.fv
#   return(r/sqrt((1+obj$sigma.fv)*obj$mu.fv))
#   #resid(obj, type = "pearson")
# 
# }
# s.fun <- function(n, obj) rNBII(n, obj$mu.fv, obj$sigma.fv)
# f.fun <- function(y) gamlss(count ~ soil_dry_mass, family = NBII(),
#                             data = spider_df)
# 
# 
# my_hnp_nblin<-hnp(fit_nb_lin, newclass = TRUE, diagfun = d.fun, simfun = s.fun,  fitfun = f.fun, resid.type = pearson, main ="negative binomial lin")

BIC(fit_p)
BIC(fit_nb)
BIC(fit_zip)
BIC(fit_zinb)
BIC(fit_nblin)



##########New code 12/02/2024 getting the estimates from hnp 


rm(list = ls(all.names = TRUE))

library(tidyverse)
library(gridExtra)
library(pscl)
#loading the data 
spider_df <- read_csv("https://raw.githubusercontent.com/rafamoral/courses/main/intro_glm/data/spider.csv")


fit_p <- glm(count ~ soil_dry_mass,
             family = poisson,
             data = spider_df)
fit_qp <- glm(count ~ soil_dry_mass,
              family = quasipoisson,
              data = spider_df)
fit_zip <- zeroinfl(count ~ soil_dry_mass,
                    dist = "poisson",
                    data = spider_df)
fit_zinb <- zeroinfl(count ~ soil_dry_mass,
                     dist = "negbin",
                     data = spider_df)
fit_nb <- glm.nb(count ~ soil_dry_mass,
                 data = spider_df)

fit_nblin <- gamlss(count ~ soil_dry_mass,
                    data = spider_df)

####getting the bic values 

bic_poisson <- BIC(fit_p)

my_distance<-function(hnp_obj) {
  r <-hnp_obj$residuals
  m<-hnp_obj$median
  


distance_p_l2 <- distance_nb_l2 <- distance_qp_l2 <- distance_zinb_l2 <- distance_zip_l2 <- distance_nblin_l2 <-  list()

distance_p_l1 <- distance_nb_l1 <- distance_qp_l1 <- distance_zinb_l1 <- distance_zip_l1 <- distance_nblin_l1 <-  list()

for(i in 1:100){
  hnp_p <- hnp(fit_p, main = "Poisson",xlab = "Half normal scores", plot = FALSE, resid.type = "pearson")
  hnp_qp <- hnp(fit_qp, main = "Quasi", xlab = "Half normal scores", plot = FALSE, resid.type = "pearson")
  hnp_zip <- hnp(fit_zip, main = "ZIP", xlab = "Half normal scores", plot = FALSE, resid.type = "pearson")
  hnp_zinb <- hnp(fit_zinb, main = "ZINB",xlab = "Half normal scores", plot = FALSE, resid.type = "pearson")
  hnp_nb <- hnp(fit_nb, main = "NB-Quad",xlab = "Half normal scores", plot = FALSE, resid.type = "pearson")
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
    
    spider_df$resp <- resp
    try_fit <- try(gamlss(resp ~ soil_dry_mass, family = NBII), silent = TRUE)
    n_try <- 1
    while(class(try_fit) == "try-error" & n_try < 20) {
      new_response <- sfun(1, object)
      try_fit <- try(gamlss(new_response ~ soil_dry_mass, family = NBII), silent = TRUE)
      n_try <- n_try + 1
    }
    return(try_fit)
    
  }
    
  
  # hnp_nblin <- hnp(fit_nblin, newclass = TRUE,diagfun = dfun, fitfun = ffun, simfun = sfun, main = "NB-lin",
                   # xlab = "Half normal scores", plot = FALSE, resid.type = "pearson")
  
  distance_p_l2[i] <- my_distance(hnp_p)[1]
  distance_nb_l2[i] <- my_distance(hnp_nb)[1]
  distance_qp_l2[i] <- my_distance(hnp_qp)[1]
  distance_zinb_l2[i] <- my_distance(hnp_zinb)[1]
  distance_zip_l2[i] <- my_distance(hnp_zip)[1]
  # distance_nblin_l2[i] <- my_distance(hnp_nblin)[1]
  
  distance_p_l1[i] <- my_distance(hnp_p)[2]
  distance_nb_l1[i] <- my_distance(hnp_nb)[2]
  distance_qp_l1[i] <- my_distance(hnp_qp)[2]
  distance_zinb_l1[i] <- my_distance(hnp_zinb)[2]
  distance_zip_l1[i] <- my_distance(hnp_zip)[2]
  # distance_nblin_l1[i] <- my_distance(hnp_nblin)[2]
  
}
distance_l2 <- tibble(distance_p_l2, distance_nb_l2, distance_qp_l2, distance_zinb_l2, distance_zip_l2)
distance_l1 <- tibble(distance_p_l1, distance_nb_l1, distance_qp_l1, distance_zinb_l1, distance_zip_l1)
# distance_l2 <- as.data.frame(distance_l2)
summary(unlist(distance_l2$distance_p_l2))
IQR(unlist(distance_l2$distance_p_l2))
sd(unlist(distance_l2$distance_p_l2))
summary(unlist(distance_l2$distance_nb_l2))
IQR(unlist(distance_l2$distance_nb_l2))
sd(unlist(distance_l2$distance_nb_l2))
summary(unlist(distance_l2$distance_qp_l2))
IQR(unlist(distance_l2$distance_qp_l2))
sd(unlist(distance_l2$distance_qp_l2))
summary(unlist(distance_l2$distance_zinb_l2))
IQR(unlist(distance_l2$distance_zinb_l2))
sd(unlist(distance_l2$distance_zinb_l2))
summary(unlist(distance_l2$distance_zip_l2))
IQR(unlist(distance_l2$distance_zip_l2))
sd(unlist(distance_l2$distance_zip_l2))

####For the L1 norm 


summary(unlist(distance_l1$distance_p_l1))
IQR(unlist(distance_l1$distance_p_l1))
sd(unlist(distance_l1$distance_p_l1))
summary(unlist(distance_l1$distance_nb_l1))
IQR(unlist(distance_l1$distance_nb_l1))
sd(unlist(distance_l1$distance_nb_l1))
summary(unlist(distance_l1$distance_qp_l1))
IQR(unlist(distance_l1$distance_qp_l1))
sd(unlist(distance_l1$distance_qp_l1))
summary(unlist(distance_l1$distance_zinb_l1))
IQR(unlist(distance_l1$distance_zinb_l1))
sd(unlist(distance_l1$distance_zinb_l1))
summary(unlist(distance_l1$distance_zip_l1))
IQR(unlist(distance_l1$distance_zip_l1))
sd(unlist(distance_l1$distance_zip_l1))



library(mvabund)
library(gamlss)
library(hnp)

data(spider)

d <- data.frame(y = spider$abund[,1],
                x = spider$x$soil.dry)

fit <- glm(y ~ x, family = quasipoisson, data = d)
hnp(fit, resid.type = "pearson")

fit2 <- gamlss(y ~ x, family = NBII, data = d)
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
    d$resp <- resp
    fit <- gamlss(resp ~ x, family = NBII, data = d)
    return(fit)
  }
  
  
  hnp_fit2 <- hnp(fit2, newclass = TRUE,
                  diagfun = dfun, fitfun = ffun, simfun = sfun)
  
  distance_fit2_l2[i] <- my_distance(hnp_fit2)[1]
  distance_fit2_l1[i] <- my_distance(hnp_fit2)[2]
  
}
summary(unlist(distance_fit2_l2))
IQR(unlist(distance_fit2_l2))
sd(unlist(distance_fit2_l2))

summary(unlist(distance_fit2_l1))
IQR(unlist(distance_fit2_l1))
sd(unlist(distance_fit2_l1))



