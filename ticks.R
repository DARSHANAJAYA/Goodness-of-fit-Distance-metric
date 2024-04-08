###Loading the libraries
library(tidyverse)
library(hnp)
library(lme4)
library(gamlss)
library(Hmisc)

#se <- function(x) sd(x)/sqrt(length(x))

###Loading the data 
tick <- read_table("Tick_field_NEPs_Marcos.txt")
### converting the variable types to factor levels 
tick$block <- as.factor(tick$block)
levels(tick$block) <- paste("Block", levels(tick$block))
tick$treat <- as.factor(tick$treat)
levels(tick$treat)[2] <- "Nematode-treated"
tick$x <- factor(1:nrow(tick))


###Exploratory plots 
png("fig_ticks.png", res = 800, units = "in", w = 6, h = 3)
tick %>%
  ggplot(aes(x = week, y = ticks, col = treat)) +
  theme_bw() +
  geom_point() +
  geom_line() +
  facet_wrap(~ block) +
  ylab("Number of ticks recovered") +
  xlab("Week") +
  labs(col = "Treatment")
dev.off()

fit_normal <- lm(ticks ~ block + factor(week) * treat,
                 data = tick)#Fitting normal model
hnp(fit_normal)#half normal plot
anova(fit_normal)#anova for thr normal model
###fiting the other. count data models and calculating the hnp objects and doing anova for the corresponding projects 
fit_poisson <- glm(ticks ~ block + factor(week) * treat,
                   family = poisson,
                   data = tick)
hnp(fit_poisson)
anova(fit_poisson, test = "Chisq")

fit_nb <- glm.nb(ticks ~ block + factor(week) * treat,
                 data = tick)
fit_nb2 <- glm.nb(ticks ~ block + factor(week) + treat,
                  data = tick)
fit_nb3 <- glm.nb(ticks ~ block + factor(week),
                  data = tick)
fit_nb4 <- glm.nb(ticks ~ block + treat,
                  data = tick)
hnp(fit_nb)
anova(fit_nb, fit_nb2, fit_nb3)
anova(fit_nb, fit_nb2, fit_nb4)

fit_poisson_mixed <- glmer(ticks ~ block + factor(week) * treat + (1 | treat : block),
                           data = tick, family = poisson)
hnp(fit_poisson_mixed)

fit_poisson_normal <- glmer(ticks ~ block + factor(week) * treat + (1 | treat : block) + (1 | x),
                            data = tick, family = poisson)
hnp(fit_poisson_normal)

drop1(fit_poisson_normal, test = "Chisq")

fit_poisson_normal2 <- glmer(ticks ~ block + factor(week) + treat + (1 | treat : block) + (1 | x),
                             data = tick, family = poisson)
drop1(fit_poisson_normal2, test = "Chisq")

## making the latextable with mean and variance
tick %>%
  group_by(treat) %>%
  summarise(mean = mean(ticks),
            se = se(ticks),
            var = var(ticks))

tick %>%
  group_by(week) %>%
  summarise(mean = mean(ticks),
            se = se(ticks),
            var = var(ticks))

tick %>%
  group_by(treat, week) %>%
  summarise(mean = mean(ticks),
            var = var(ticks)) %>%
  latexTabular %>%
  cat

png("ticks_normal_resid.png", res = 800, units = "in", w = 6, h = 4)
tibble(residuals = residuals(fit_normal),
       fitted = fitted(fit_normal)) %>%
  ggplot(aes(x = fitted, y = residuals)) +
  theme_bw() +
  geom_abline(intercept = 0, slope = 0, lty = 2) +
  geom_point() +
  xlab("Fitted values") +
  ylab("Raw residuals")
dev.off()
####hnp model objects and half normal plots 

set.seed(2023)
hnp1 <- hnp(fit_normal)
hnp2 <- hnp(fit_poisson)
hnp3 <- hnp(fit_nb)

png("ticks_hnp.png", res = 800, units = "in", w = 9, h = 3)
par(mfrow = c(1,3))
plot(hnp1, xlab = "Half-normal scores", ylab = "Raw residuals", main = "(a)")
plot(hnp2, xlab = "Half-normal scores", ylab = "Deviance residuals", main = "(b)")
plot(hnp3, xlab = "Half-normal scores", ylab = "Deviance residuals", main = "(c)")
dev.off()


