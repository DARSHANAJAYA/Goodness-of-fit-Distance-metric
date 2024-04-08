rm(list = ls(all.names = TRUE)) #Clear the workspace 

###Loading the libraries 
library(readxl)
library(tidyr)
library(ggplot2)
library(fitdistrplus)
library(actuar)
library(hnp)

###Loading the data 

plant <- read_excel("All plant traits and nematode per gram dry weight Heterodera.xlsx")

###Exploratory plot
p<-ggplot(data = plant, aes(x = dry_weight, color = treat)) +
  xlim(-0.6, 0.8)+
  geom_density()+
  labs(color = "treatment")+
  theme_classic()


###Fitting the models and calculating the haf normal objects for the corresponding model

fit_normal <- lm(dry_weight~treat, data = plant)
hnp_normal <- hnp(fit_normal, resid.type ="pearson", xlab = "Half normal scores", ylab ="Pearson residuals",main ="(a) Normal model")
#fit <- fitdist(plant$dry_weight, distr = "gamma", method = "mle")
fit_gamma <- glm(dry_weight~treat,data = plant, family = Gamma)
hnp_gamma <- hnp(fit_gamma, resid.type ="pearson", xlab = "Half normal scores", ylab ="Pearson residuals",main ="(b) Gamma model")

fit_invgauss <- glm(dry_weight~treat,data = plant, family = inverse.gaussian)
hnp_invgauss <- hnp(fit_invgauss, resid.type ="pearson", xlab = "Half normal scores", ylab ="Pearson residuals",main ="(c) Inverse Gaussian model")

 
par(mfrow = c(1, 3))
###half normal plots 

plot(hnp_normal,xlab = "Half normal scores", ylab ="Pearson residuals",main ="(a) Normal model",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
plot(hnp_gamma,xlab = "Half normal scores", ylab ="Pearson residuals",main ="(b) Gamma model",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
plot(hnp_invgauss,xlab = "Half normal scores", ylab ="Pearson residuals",main ="(c) Inverse Gaussian model",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)

