####Code for poisson 
rm(list = ls(all.names = TRUE))
set.seed(123)
library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)
library(ggpubr)
theme_set(theme_pubr())

L2 <- as.data.frame(rep("p = 2", times = 44865))
colnames(L2)<- c("L")
L1 <- as.data.frame(rep("p = 1", times = 44865))
colnames(L1)<- c("L")
L_norm <- as.data.frame(rbind(L2, L1))

twenty <- rep(20, times =997)
fifty <- rep(50, times =997)
hundred <- rep(100, times = 997)
sample_size <- rep(c(twenty,fifty,hundred), times = 30)
sample_size <- as.data.frame(sample_size)


######This for the frequency of distances 

datafinalpois13 <- read.csv("datafinalupdatedpois13.csv")
datafinal1 <- datafinalpois13 %>% 
  slice(rep(1:n(), 15))
datafinalpois4648 <-  read.csv("poisson/datafinalupdatedpois4648.csv")
datafinal2 <- datafinalpois4648 %>% 
  slice(rep(1:n(), 15))

datafinalpois <- rbind(datafinal1, datafinal2)
model1 <-  rep("Poisson", times = nrow(datafinalpois))

#run this only after running the Lnorm and the sample size 

datafinalpois1 <- data.frame(cbind(base_value = datafinalpois$base_value_poisson, model1, L_norm, sample_size)) ### this variable is for the overall distance




datafinalnegbin13 <- read.csv("poisson/datafinalupdatednbquad13.csv")
datafinal1 <- datafinalnegbin13 %>% 
  slice(rep(1:n(), 15))
datafinalnegbin4648 <-  read.csv("poisson/datafinalupdatednbquad4648.csv")
datafinal2 <- datafinalnegbin4648 %>% 
  slice(rep(1:n(), 15))

datafinalnegquad <- rbind(datafinal1, datafinal2)
model1 <-  rep("NB-quad", times = nrow(datafinalnegquad))
datafinalnegquad1 <- data.frame(cbind(base_value = datafinalnegquad$base_value_negbinquad, model1, L_norm, sample_size))

datafinalnblin13 <- read.csv("poisson/datafinalupdatednblin13.csv")
datafinal1 <- datafinalnblin13 %>% 
  slice(rep(1:n(), 15))
datafinalnblin4648 <-  read.csv("/Share/localshare/djayakumari/poisson/datafinalupdatednblin4648.csv")
datafinal2 <- datafinalnblin4648 %>% 
  slice(rep(1:n(), 15))

datafinalnblin <- rbind(datafinal1, datafinal2)

model1 <-  rep("NB-lin", times = nrow(datafinalnblin))
datafinalnblin1 <- data.frame(cbind(base_value = datafinalnblin$base_value_negbinlin, model1, L_norm, sample_size))



datafinalquasi13 <- read.csv("poisson/datafinalupdatedquasi13.csv")
datafinal1 <- datafinalquasi13 %>% 
  slice(rep(1:n(), 15))
datafinalquasi4648 <-  read.csv("poisson/datafinalupdatedquasi4648.csv")
datafinal2 <- datafinalquasi4648 %>% 
  slice(rep(1:n(), 15))

datafinalquasi <- rbind(datafinal1, datafinal2)
model1 <-  rep("Quasi", times = nrow(datafinalquasi))
datafinalquasi1 <- data.frame(cbind(base_value = datafinalquasi$base_value_quasi, model1, L_norm, sample_size))

# datafinalzeroinflnegbin13 <- read.csv("/Users/darshanaj/Documents/individual project/sim_new/dataframes/zeroinflnegbin/datafinalupdatedzeroinflnb13.csv")
# datafinal1 <- datafinalzeroinflnegbin13 %>% 
#   slice(rep(1:n(), 15))
# datafinalzeroinflnegbin4648 <-  read.csv("/Users/darshanaj/Documents/individual project/sim_new/dataframes/zeroinflnegbin/datafinalupdatedzeroinflnb4648.csv")
# datafinal2 <- datafinalzeroinflnegbin4648 %>% 
#   slice(rep(1:n(), 15))
# 
# datafinalzinbh <- rbind(datafinal1, datafinal2)


base_distance <- as.data.frame(cbind(datafinalpois$base_value_poisson, datafinalquasi$base_value_quasi,
                                     datafinalnblin$base_value_negbinlin,datafinalnegquad$base_value_negbinquad))


names(base_distance)  <- c("Poisson", "Quasi", "NB-lin", "NB-quad")




base_distance <- cbind(base_distance, L_norm, sample_size)
base_distance$sample_size <- as.factor(base_distance$sample_size)
base_distance$L <- as.factor(base_distance$L)

base_distance <- na.omit(base_distance)
base_distance$min_model <- names(base_distance[,1:4])[apply(base_distance[,1:4], MARGIN = 1, FUN = which.min)]

base_distance$min_model <-  factor(base_distance$min_model,     # Reorder factor levels
                                   c("Poisson", "Quasi", "NB-lin", "NB-quad"))
summary = base_distance %>% group_by(L, sample_size,min_model) %>%
  tally %>%
  group_by(L, sample_size) %>%
  mutate(pct = n/sum(n))

p6 <- ggplot(summary, aes(min_model)) +
  geom_bar(aes(y = pct, fill = factor(min_model)), stat = "identity", show.legend = FALSE) +
  theme_pubclean() +
  facet_wrap(~ sample_size + L, nrow = 4, labeller =
               labeller(
                 sample_size = ~ paste("sample size:", .),
                 L = ~ paste(.),
                 .multi_line = FALSE)) +
  labs(x = "Model", y = "% smallest d") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = c("Poisson" = "#1B9E77", "Quasi" = "#D95F02",
                               "NB-lin" = "#7570B3", "NB-quad" = "#E7298A"))

p6


#### For the overall distance plot


overall_dist <- rbind(datafinalpois1, datafinalquasi1, datafinalnblin1, datafinalnegquad1)
overall_dist$sample_size <- as.factor(overall_dist$sample_size)
overall_dist$L <- as.factor(overall_dist$L)
overall_dist <- na.omit(overall_dist)

overall_dist$model1 <-  factor(overall_dist$model1,     # Reorder factor levels
                               c("Poisson", "Quasi", "NB-lin", "NB-quad"))

p1 <- overall_dist %>% ggplot(aes(x = sample_size, y = log(base_value))) +
  theme_bw() +
  geom_boxplot(aes(color = model1)) +
  facet_wrap(~L)+ 
  labs(x = "Sample size", y = "ln(d)")+
  ggtitle(expression("Parent model: Poisson"))+
  theme(legend.position = "top")+
  scale_color_manual( name = "Fitted models",values=c("Poisson" = "#1B9E77", "Quasi" ="#D95F02", 
                                                      "NB-lin"= "#7570B3", "NB-quad" = "#E7298A"))
p1
####for the bic 

bic_poisson <- data.frame(cbind(bic_poisson = datafinalpois$bic_poisson,
                                bic_negbinlin = datafinalnblin$bic_nblin, bic_negbinquad = datafinalnegquad$bic_nbquad, L_norm, sample_size))
names(bic_poisson) <- c("Poisson","NB-lin", "NB-quad", "L_norm","sample_size" )

bic_poisson$sample_size <- as.factor(bic_poisson$sample_size)
#bic_poisson$L <- as.factor(bic_poisson$L)

bic_poisson <- na.omit(bic_poisson)
bic_poisson$min_model <- names(bic_poisson[,1:3])[apply(bic_poisson[,1:3], MARGIN = 1, FUN = which.min)]
bic_poisson$min_model <- as.factor(bic_poisson$min_model)
p7 <-ggplot(bic_poisson, aes(fct_infreq(min_model)))+
  geom_bar(aes(y = (..count..)/sum(..count..), fill =min_model), show.legend = FALSE) +
  theme_pubclean()+
  labs(x = "Model", y = "% smallest BIC")+ 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))+
  
  scale_fill_manual(values=c("Poisson" = "#1B9E77", "Quasi" ="#D95F02", 
                             "NB-lin"= "#7570B3", "NB-quad" = "#E7298A"))+
  coord_fixed(ratio =1.7)


p7
ggarrange(p1,                                                 # First row with scatter plot
          ggarrange(p6, p7, nrow = 2, labels = c("b", "c")), # Second row with box and dot plots
          ncol = 2, 
          labels = "a"                                        # Labels of the scatter plot
) 

###code for negbin
rm(list = ls(all.names = TRUE))
set.seed(123)
library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)
library(ggpubr)
theme_set(theme_pubr())

L2 <- as.data.frame(rep("p = 2", times = 44910))
colnames(L2)<- c("L")
L1 <- as.data.frame(rep("p = 1", times = 44910))
colnames(L1)<- c("L")
L_norm <- as.data.frame(rbind(L2, L1))

twenty <- rep(20, times =998)
fifty <- rep(50, times =998)
hundred <- rep(100, times = 998)
sample_size <- rep(c(twenty,fifty,hundred), times = 30)
sample_size <- as.data.frame(sample_size)

######This for the frequency of distances 

datafinalpois13 <- read.csv("dataframenegbinm/datafinalupdatedpois13.csv")
datafinal1 <- datafinalpois13 %>% 
  slice(rep(1:n(), 15))
datafinalpois4648 <-  read.csv("dataframenegbinm/datafinalupdatedpois4648.csv")
datafinal2 <- datafinalpois4648 %>% 
  slice(rep(1:n(), 15))

datafinalpois <- rbind(datafinal1, datafinal2)
model1 <-  rep("Poisson", times = nrow(datafinalpois))

#run this only after running the Lnorm and the sample size 

datafinalpois1 <- data.frame(cbind(base_value = datafinalpois$base_value_poisson, model1, L_norm, sample_size)) ### this variable is for the overall distance




datafinalnegbin13 <- read.csv("dataframenegbinm/datafinalupdatednbquad13.csv")
datafinal1 <- datafinalnegbin13 %>% 
  slice(rep(1:n(), 15))
datafinalnegbin4648 <-  read.csv("dataframenegbinm/datafinalupdatednbquad4648.csv")
datafinal2 <- datafinalnegbin4648 %>% 
  slice(rep(1:n(), 15))

datafinalnegquad <- rbind(datafinal1, datafinal2)
model1 <-  rep("NB-quad", times = nrow(datafinalnegquad))
datafinalnegquad1 <- data.frame(cbind(base_value = datafinalnegquad$base_value_negbinquad, model1, L_norm, sample_size))

datafinalnblin13 <- read.csv("dataframenegbinm/datafinalupdatednblin13.csv")
datafinal1 <- datafinalnblin13 %>% 
  slice(rep(1:n(), 15))
datafinalnblin4648 <-  read.csv("dataframenegbinm/datafinalupdatednblin4648.csv")
datafinal2 <- datafinalnblin4648 %>% 
  slice(rep(1:n(), 15))

datafinalnblin <- rbind(datafinal1, datafinal2)

model1 <-  rep("NB-lin", times = nrow(datafinalnblin))
datafinalnblin1 <- data.frame(cbind(base_value = datafinalnblin$base_value_negbinlin, model1, L_norm, sample_size))



datafinalquasi13 <- read.csv("dataframenegbinm/datafinalupdatedquasi13.csv")
datafinal1 <- datafinalquasi13 %>% 
  slice(rep(1:n(), 15))
datafinalquasi4648 <-  read.csv("new/dataframenegbinm/datafinalupdatedquasi4648.csv")
datafinal2 <- datafinalquasi4648 %>% 
  slice(rep(1:n(), 15))

datafinalquasi <- rbind(datafinal1, datafinal2)
model1 <-  rep("Quasi", times = nrow(datafinalquasi))
datafinalquasi1 <- data.frame(cbind(base_value = datafinalquasi$base_value_quasi, model1, L_norm, sample_size))

# datafinalzeroinflnegbin13 <- read.csv("/Users/darshanaj/Documents/individual project/sim_new/dataframes/zeroinflnegbin/datafinalupdatedzeroinflnb13.csv")
# datafinal1 <- datafinalzeroinflnegbin13 %>% 
#   slice(rep(1:n(), 15))
# datafinalzeroinflnegbin4648 <-  read.csv("/Users/darshanaj/Documents/individual project/sim_new/dataframes/zeroinflnegbin/datafinalupdatedzeroinflnb4648.csv")
# datafinal2 <- datafinalzeroinflnegbin4648 %>% 
#   slice(rep(1:n(), 15))
# 
# datafinalzinbh <- rbind(datafinal1, datafinal2)

base_distance <- as.data.frame(cbind(datafinalpois$base_value_poisson, datafinalquasi$base_value_quasi,
                                     datafinalnblin$base_value_negbinlin,datafinalnegquad$base_value_negbinquad))



names(base_distance)  <- c("Poisson", "Quasi", "NB-lin", "NB-quad")





base_distance <- cbind(base_distance, L_norm, sample_size)
base_distance$sample_size <- as.factor(base_distance$sample_size)
base_distance$L <- as.factor(base_distance$L)

base_distance <- na.omit(base_distance)
base_distance$min_model <- names(base_distance[,1:4])[apply(base_distance[,1:4], MARGIN = 1, FUN = which.min)]

base_distance$min_model <-  factor(base_distance$min_model,     # Reorder factor levels
                                   c("Poisson", "Quasi", "NB-lin", "NB-quad"))

summary = base_distance %>% group_by(L, sample_size,min_model) %>%
  tally %>%
  group_by(L, sample_size) %>%
  mutate(pct = n/sum(n))

p6 <- ggplot(summary, aes(min_model)) +
  geom_bar(aes(y = pct, fill = factor(min_model)), stat = "identity", show.legend = FALSE) +
  theme_pubclean() +
  facet_wrap(~ sample_size + L, nrow = 4, labeller =
               labeller(
                 sample_size = ~ paste("sample size:", .),
                 L = ~ paste(.),
                 .multi_line = FALSE)) +
  labs(x = "Model", y = "% smallest d") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = c("Poisson" = "#1B9E77", "Quasi" = "#D95F02",
                               "NB-lin" = "#7570B3", "NB-quad" = "#E7298A"))

p6
#### For the overall distance plot


overall_dist <- rbind(datafinalpois1, datafinalquasi1, datafinalnblin1, datafinalnegquad1)
overall_dist$sample_size <- as.factor(overall_dist$sample_size)
overall_dist$L <- as.factor(overall_dist$L)
overall_dist <- na.omit(overall_dist)
overall_dist$model1 <- as.factor(overall_dist$model1)
overall_dist$model1 <-  factor(overall_dist$model1,     # Reorder factor levels
                               c("Poisson", "Quasi", "NB-lin", "NB-quad"))

p1 <- overall_dist %>% ggplot(aes(x = sample_size, y = log(base_value))) +
  theme_bw() +
  geom_boxplot(aes(color = model1)) +
  facet_wrap(~L)+ 
  labs(x = "Sample size", y = "ln(d)")+
  ggtitle(expression("Parent model: NB-quad("*zeta~ "= 0.5)"))+
  theme(legend.position = "top")+
  scale_color_manual( name = "Fitted models",values=c("Poisson" = "#1B9E77", "Quasi" ="#D95F02", 
                                                      "NB-lin"= "#7570B3", "NB-quad" = "#E7298A"))
p1
bic_negbinm <- data.frame(cbind(bic_poisson = datafinalpois$bic_poisson,
                                bic_negbinlin = datafinalnblin$bic_nblin, bic_negbinquad = datafinalnegquad$bic_nbquad, L_norm, sample_size))
names(bic_negbinm) <- c("Poisson","NB-lin", "NB-quad", "L_norm","sample_size" )
bic_negbinm$sample_size <- as.factor(bic_negbinm$sample_size)
bic_negbinm$L <- as.factor(bic_negbinm$L)

bic_negbinm <- na.omit(bic_negbinm)
bic_negbinm$min_model <- names(bic_negbinm[,1:3])[apply(bic_negbinm[,1:3], MARGIN = 1, FUN = which.min)]
bic_negbinm$min_model<-  factor(bic_negbinm$min_model,     # Reorder factor levels
                                c("Poisson", "NB-lin", "NB-quad"))
p7 <-ggplot(bic_negbinm, aes(min_model)) +
  geom_bar(aes(y = (..count..)/sum(..count..), fill =min_model), show.legend = FALSE) +
  theme_pubclean()+ 
  labs(x = "Model", y = "% smallest BIC")+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))+
  scale_fill_manual(values=c("Poisson" = "#1B9E77", "Quasi" ="#D95F02", 
                             "NB-lin"= "#7570B3", "NB-quad" = "#E7298A"))+
  coord_fixed(ratio =2.2)
p7
# finalplot <- ggarrange(p1, p6, p7, ncol = 3, nrow = 1)
# annotate_figure(finalplot, top = text_grob("Parent model:Negbin-quad(strong overdispersion)", 
#                                            color = "red", face = "bold", size = 14))
ggarrange(p1,                                                 # First row with scatter plot
          ggarrange(p6, p7, nrow = 2, labels = c("b", "c")), # Second row with box and dot plots
          ncol = 2, 
          labels = "a"                                        # Labels of the scatter plot
) 

### code for negbinlin

rm(list = ls(all.names = TRUE))
set.seed(123)
library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)
library(ggpubr)
theme_set(theme_pubr())

######This for the frequency of distances 

L2 <- as.data.frame(rep("p = 2", times = 44910))
colnames(L2)<- c("L")
L1 <- as.data.frame(rep("p = 1", times = 44910))
colnames(L1)<- c("L")
L_norm <- as.data.frame(rbind(L2, L1))

twenty <- rep(20, times =998)
fifty <- rep(50, times =998)
hundred <- rep(100, times = 998)
sample_size <- rep(c(twenty,fifty,hundred), times = 30)
sample_size <- as.data.frame(sample_size)

datafinalpois13 <- read.csv("dataframenegbinhlin/datafinalupdatedpois13.csv")
datafinal1 <- datafinalpois13 %>% 
  slice(rep(1:n(), 15))
datafinalpois4648 <-  read.csv("dataframenegbinhlin/datafinalupdatedpois4648.csv")
datafinal2 <- datafinalpois4648 %>% 
  slice(rep(1:n(), 15))

datafinalpois <- rbind(datafinal1, datafinal2)
model1 <-  rep("Poisson", times = nrow(datafinalpois))

#run this only after running the Lnorm and the sample size 

datafinalpois1 <- data.frame(cbind(base_value = datafinalpois$base_value_poisson, model1, L_norm, sample_size)) ### this variable is for the overall distance




datafinalnegbin13 <- read.csv("dataframenegbinhlin/datafinalupdatednbquad13.csv")
datafinal1 <- datafinalnegbin13 %>% 
  slice(rep(1:n(), 15))
datafinalnegbin4648 <-  read.csv("dataframenegbinhlin/datafinalupdatednbquad4648.csv")
datafinal2 <- datafinalnegbin4648 %>% 
  slice(rep(1:n(), 15))

datafinalnegquad <- rbind(datafinal1, datafinal2)
model1 <-  rep("NB-quad", times = nrow(datafinalnegquad))
datafinalnegquad1 <- data.frame(cbind(base_value = datafinalnegquad$base_value_negbinquad, model1, L_norm, sample_size))

datafinalnblin13 <- read.csv("dataframenegbinhlin/datafinalupdatednblin13.csv")
datafinal1 <- datafinalnblin13 %>% 
  slice(rep(1:n(), 15))
datafinalnblin4648 <-  read.csv("dataframenegbinhlin/datafinalupdatednblin4648.csv")
datafinal2 <- datafinalnblin4648 %>% 
  slice(rep(1:n(), 15))

datafinalnblin <- rbind(datafinal1, datafinal2)

model1 <-  rep("NB-lin", times = nrow(datafinalnblin))
datafinalnblin1 <- data.frame(cbind(base_value = datafinalnblin$base_value_negbinlin, model1, L_norm, sample_size))



datafinalquasi13 <- read.csv("dataframenegbinhlin/datafinalupdatedquasi13.csv")
datafinal1 <- datafinalquasi13 %>% 
  slice(rep(1:n(), 15))
datafinalquasi4648 <-  read.csv("dataframenegbinhlin/datafinalupdatedquasi4648.csv")
datafinal2 <- datafinalquasi4648 %>% 
  slice(rep(1:n(), 15))

datafinalquasi <- rbind(datafinal1, datafinal2)
model1 <-  rep("Quasi", times = nrow(datafinalquasi))
datafinalquasi1 <- data.frame(cbind(base_value = datafinalquasi$base_value_quasi, model1, L_norm, sample_size))

# datafinalzeroinflnegbin13 <- read.csv("/Users/darshanaj/Documents/individual project/sim_new/dataframes/zeroinflnegbin/datafinalupdatedzeroinflnb13.csv")
# datafinal1 <- datafinalzeroinflnegbin13 %>% 
#   slice(rep(1:n(), 15))
# datafinalzeroinflnegbin4648 <-  read.csv("/Users/darshanaj/Documents/individual project/sim_new/dataframes/zeroinflnegbin/datafinalupdatedzeroinflnb4648.csv")
# datafinal2 <- datafinalzeroinflnegbin4648 %>% 
#   slice(rep(1:n(), 15))
# 
# datafinalzinbh <- rbind(datafinal1, datafinal2)

base_distance <- as.data.frame(cbind(datafinalpois$base_value_poisson, datafinalquasi$base_value_quasi,
                                     datafinalnblin$base_value_negbinlin,datafinalnegquad$base_value_negbinquad))


names(base_distance)  <- c("Poisson", "Quasi", "NB-lin", "NB-quad")


L2 <- as.data.frame(rep("p = 2", times = 44910))
colnames(L2)<- c("L")
L1 <- as.data.frame(rep("p = 1", times = 44910))
colnames(L1)<- c("L")
L_norm <- as.data.frame(rbind(L2, L1))

twenty <- rep(20, times =998)
fifty <- rep(50, times =998)
hundred <- rep(100, times = 998)
sample_size <- rep(c(twenty,fifty,hundred), times = 30)
sample_size <- as.data.frame(sample_size)

base_distance <- cbind(base_distance, L_norm, sample_size)
base_distance$sample_size <- as.factor(base_distance$sample_size)
base_distance$L <- as.factor(base_distance$L)

base_distance <- na.omit(base_distance)
base_distance$min_model <- names(base_distance[,1:4])[apply(base_distance[,1:4], MARGIN = 1, FUN = which.min)]

base_distance$min_model <-  factor(base_distance$min_model,     # Reorder factor levels
                                   c("Poisson", "Quasi", "NB-lin", "NB-quad"))
summary = base_distance %>% group_by(L, sample_size,min_model) %>%
  tally %>%
  group_by(L, sample_size) %>%
  mutate(pct = n/sum(n))

p6 <- ggplot(summary, aes(min_model)) +
  geom_bar(aes(y = pct, fill = factor(min_model)), stat = "identity", show.legend = FALSE) +
  theme_pubclean() +
  facet_wrap(~ sample_size + L, nrow = 4, labeller =
               labeller(
                 sample_size = ~ paste("sample size:", .),
                 L = ~ paste(.),
                 .multi_line = FALSE)) +
  labs(x = "Model", y = "% smallest d") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = c("Poisson" = "#1B9E77", "Quasi" = "#D95F02",
                               "NB-lin" = "#7570B3", "NB-quad" = "#E7298A"))

p6
#### For the overall distance plot


overall_dist <- rbind(datafinalpois1, datafinalquasi1, datafinalnblin1, datafinalnegquad1)
overall_dist$sample_size <- as.factor(overall_dist$sample_size)
overall_dist$L <- as.factor(overall_dist$L)
overall_dist <- na.omit(overall_dist)
overall_dist$model1 <-  factor(overall_dist$model1,     # Reorder factor levels
                               c("Poisson", "Quasi", "NB-lin", "NB-quad"))

p1 <- overall_dist %>% ggplot(aes(x = sample_size, y = log(base_value))) +
  theme_bw() +
  geom_boxplot(aes(color = model1)) +
  facet_wrap(~L)+ 
  labs(x = "Sample size", y = "ln(d)")+
  ggtitle(expression("Parent model: NB-lin("*phi~ "= 7)"))+
  theme(legend.position = "top")+
  scale_color_manual( name = "Fitted models",values=c("Poisson" = "#1B9E77", "Quasi" ="#D95F02", 
                                                      "NB-lin"= "#7570B3", "NB-quad" = "#E7298A"))
p1
####for the bic 
bic_negbinhlin <- data.frame(cbind(bic_poisson = datafinalpois$bic_poisson,
                                   bic_negbinlin = datafinalnblin$bic_nblin, bic_negbinquad = datafinalnegquad$bic_nbquad, L_norm, sample_size))

names(bic_negbinhlin) <- c("Poisson","NB-lin", "NB-quad", "L_norm","sample_size" )
bic_negbinhlin$sample_size <- as.factor(bic_negbinhlin$sample_size)
bic_negbinhlin$L <- as.factor(bic_negbinhlin$L)

bic_negbinhlin <- na.omit(bic_negbinhlin)
bic_negbinhlin$min_model <- names(bic_negbinhlin[,1:3])[apply(bic_negbinhlin[,1:3], MARGIN = 1, FUN = which.min)]
bic_negbinhlin$min_model<-  factor(bic_negbinhlin$min_model,     # Reorder factor levels
                                   c("Poisson", "NB-lin", "NB-quad"))
p7 <-ggplot(bic_negbinhlin, aes(min_model)) +
  geom_bar(aes(y = (..count..)/sum(..count..), fill =min_model), show.legend = FALSE) +
  theme_pubclean()+
  labs(x = "Model", y = "% smallest BIC")+ 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))+
  scale_fill_manual(values=c("Poisson" = "#1B9E77", "Quasi" ="#D95F02", 
                             "NB-lin"= "#7570B3", "NB-quad" = "#E7298A"))+
  coord_fixed(ratio =1.5)
p7
# finalplot <- ggarrange(p1, p6, p7, ncol = 3, nrow = 1)
# annotate_figure(finalplot, top = text_grob("Parent model:Negbin-quad(strong overdispersion)", 
#                                            color = "red", face = "bold", size = 14))
ggarrange(p1,                                                 # First row with scatter plot
          ggarrange(p6, p7, nrow = 2, labels = c("b", "c")), # Second row with box and dot plots
          ncol = 2, 
          labels = "a"                                        # Labels of the scatter plot
) 


###code for zeroinflated
rm(list = ls(all.names = TRUE))
set.seed(123)
library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)
library(ggpubr)
theme_set(theme_pubr())

######This for the frequency of distances 

L2 <- as.data.frame(rep("p = 2", times = 44775))
colnames(L2)<- c("L")
L1 <- as.data.frame(rep("p = 1", times = 44775))
colnames(L1)<- c("L")
L_norm <- as.data.frame(rbind(L2, L1))

twenty <- rep(20, times =995)
fifty <- rep(50, times =995)
hundred <- rep(100, times = 995)
sample_size <- rep(c(twenty,fifty,hundred), times = 30)
sample_size <- as.data.frame(sample_size)


datafinalpois13 <- read.csv("zinph/datafinalupdatedpois13.csv")
datafinal1 <- datafinalpois13 %>% 
  slice(rep(1:n(), 15))
datafinalpois4648 <-  read.csv("zinph/datafinalupdatedpois4648.csv")
datafinal2 <- datafinalpois4648 %>% 
  slice(rep(1:n(), 15))


datafinalpois <- rbind(datafinal1, datafinal2)
model1 <-  rep("Poisson", times = nrow(datafinalpois))

#run this only after running the Lnorm and the sample size 

datafinalpois1 <- data.frame(cbind(base_value = datafinalpois$base_value_poisson, model1, L_norm, sample_size)) ### this variable is for the overall distance




datafinalnegbin13 <- read.csv("zinph/datafinalupdatednbquad13.csv")
datafinal1 <- datafinalnegbin13 %>% 
  slice(rep(1:n(), 15))
datafinalnegbin4648 <-  read.csv("zinph/datafinalupdatednbquad4648.csv")
datafinal2 <- datafinalnegbin4648 %>% 
  slice(rep(1:n(), 15))

datafinalnegquad <- rbind(datafinal1, datafinal2)
model1 <-  rep("NB-quad", times = nrow(datafinalnegquad))
datafinalnegquad1 <- data.frame(cbind(base_value = datafinalnegquad$base_value_negbinquad, model1, L_norm, sample_size))

datafinalnblin13 <- read.csv("zinph/datafinalupdatednblin13.csv")
datafinal1 <- datafinalnblin13 %>% 
  slice(rep(1:n(), 15))
datafinalnblin4648 <-  read.csv("zinph/datafinalupdatednblin4648.csv")
datafinal2 <- datafinalnblin4648 %>% 
  slice(rep(1:n(), 15))

datafinalnblin <- rbind(datafinal1, datafinal2)

model1 <-  rep("NB-lin", times = nrow(datafinalnblin))
datafinalnblin1 <- data.frame(cbind(base_value = datafinalnblin$base_value_negbinlin, model1, L_norm, sample_size))



datafinalquasi13 <- read.csv("zinph/datafinalupdatedquasi13.csv")
datafinal1 <- datafinalquasi13 %>% 
  slice(rep(1:n(), 15))
datafinalquasi4648 <-  read.csv("datafinalupdatedquasi4648.csv")
datafinal2 <- datafinalquasi4648 %>% 
  slice(rep(1:n(), 15))

datafinalquasi <- rbind(datafinal1, datafinal2)
model1 <-  rep("Quasi", times = nrow(datafinalquasi))
datafinalquasi1 <- data.frame(cbind(base_value = datafinalquasi$base_value_quasi, model1, L_norm, sample_size))

datafinalzeroinflnegbin13 <- read.csv("zinph/datafinalupdatedzeroinflnb13.csv")
datafinal1 <- datafinalzeroinflnegbin13 %>%
  slice(rep(1:n(), 15))
datafinalzeroinflnegbin4648 <-  read.csv("zinph/datafinalupdatedzeroinflnb4648.csv")
datafinal2 <- datafinalzeroinflnegbin4648 %>%
  slice(rep(1:n(), 15))

datafinalzinbh <- rbind(datafinal1, datafinal2)

model1 <-  rep("ZINB", times = nrow(datafinalzinbh))
datafinalzinbh1 <- data.frame(cbind(base_value = datafinalzinbh$base_value_zinbh, model1, L_norm, sample_size))


datafinalzeroinflpoisson13 <- read.csv("zinph/datafinalupdatedzeroinflp13.csv")
datafinal1 <- datafinalzeroinflpoisson13 %>%
  slice(rep(1:n(), 15))
datafinalzeroinflpoisson4648 <-  read.csv("zinph/datafinalupdatedzeroinflp4648.csv")
datafinal2 <- datafinalzeroinflpoisson4648 %>%
  slice(rep(1:n(), 15))

datafinalzip <- rbind(datafinal1, datafinal2)

model1 <-  rep("ZIP", times = nrow(datafinalzip))
datafinalzip1 <- data.frame(cbind(base_value = datafinalzip$base_value_zinp, model1, L_norm, sample_size))


base_distance <- as.data.frame(cbind(poisson = datafinalpois$base_value_poisson, quasi = datafinalquasi$base_value_quasi,
                                     negbinlin = datafinalnblin$base_value_negbinlin, negbinquad = datafinalnegquad$base_value_negbinquad,
                                     zip = datafinalzip$base_value_zinp,zinb = datafinalzinbh$base_value_zinbh))

names(base_distance)  <- c("Poisson", "Quasi", "NB-lin", "NB-quad", "ZIP", "ZINB")
# names(base_distance)  <- c("poisson", "negbinquad", "negbinlin", "quasi")


# names(base_distance)  <- c("poisson", "negbinquad", "negbinlin", "quasi")



base_distance <- cbind(base_distance, L_norm, sample_size)
base_distance$sample_size <- as.factor(base_distance$sample_size)
base_distance$L <- as.factor(base_distance$L)

base_distance <- na.omit(base_distance)
base_distance$min_model <- names(base_distance[,1:6])[apply(base_distance[,1:6], MARGIN = 1, FUN = which.min)]




base_distance$min_model <-  factor(base_distance$min_model,     # Reorder factor levels
                                   c("Poisson", "Quasi", "NB-lin", "NB-quad", "ZIP", "ZINB"))

summary = base_distance %>% group_by(L, sample_size,min_model) %>%
  tally %>%
  group_by(L, sample_size) %>%
  mutate(pct = n/sum(n))

p6 <- ggplot(summary, aes(min_model)) +
  geom_bar(aes(y = pct, fill = factor(min_model)), stat = "identity", show.legend = FALSE) +
  theme_pubclean() +
  facet_wrap(~ sample_size + L, nrow = 4, labeller =
               labeller(
                 sample_size = ~ paste("sample size:", .),
                 L = ~ paste(.),
                 .multi_line = FALSE)) +
  labs(x = "Model", y = "% smallest d") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values=c("Poisson" = "#1B9E77", "Quasi" ="#D95F02", 
                             "NB-lin"= "#7570B3", "NB-quad" = "#E7298A", "ZIP" = "#084594", "ZINB" = "#A6761D"))

p6

# names(base_distance)  <- c("poisson", "negbinquad", "negbinlin", "quasi")
#####for the overall distance plot


overall_dist <- rbind(datafinalpois1, datafinalquasi1, datafinalnblin1, datafinalnegquad1, datafinalzip1, datafinalzinbh1)
overall_dist$sample_size <- as.factor(overall_dist$sample_size)
overall_dist$L <- as.factor(overall_dist$L)
overall_dist <- na.omit(overall_dist)

overall_dist$model1 <-  factor(overall_dist$model1,     # Reorder factor levels
                               c("Poisson", "Quasi", "NB-lin", "NB-quad","ZIP", "ZINB"))

p1 <- overall_dist %>% ggplot(aes(x = sample_size, y = log(base_value))) +
  theme_bw() +
  geom_boxplot(aes(color = model1)) +
  facet_wrap(~L)+ 
  labs(x = "Sample size", y = "ln(d)")+
  ggtitle(expression("Parent model: ZIP("*nu~ "= 0.6)"))+
  theme(legend.position = "top")+
  scale_color_manual( name = "Fitted models",values=c("Poisson" = "#1B9E77", "Quasi" ="#D95F02", 
                                                      "NB-lin"= "#7570B3", "NB-quad" = "#E7298A", 
                                                      "ZIP" = "#084594", "ZINB" = "#A6761D" ))
p1
####for the bic 

bic_zinph <- data.frame(cbind(bic_poisson = datafinalpois$bic_poisson,
                              bic_neglin = datafinalnblin$bic_neg_lin,
                              bic_negquad = datafinalnegquad$bic_neg_quad,
                              bic_zip = datafinalzip$bic_zinp,
                              bic_zinb = datafinalzinbh$bic_zinbh,
                              L_norm, sample_size))
names(bic_zinph) <- c("Poisson", "NB-lin", "NB-quad","ZIP", "ZINB", "L_norm","sample_size" )

bic_zinph$sample_size <- as.factor(bic_zinph$sample_size)
bic_zinph$L <- as.factor(bic_zinph$L)

bic_zinph <- na.omit(bic_zinph)
bic_zinph$min_model <- names(bic_zinph[,1:5])[apply(bic_zinph[,1:5], MARGIN = 1, FUN = which.min)]
bic_zinph$min_model<-  factor(bic_zinph$min_model,     # Reorder factor levels
                              c("Poisson", "NB-lin", "NB-quad","ZIP", "ZINB"))
p7 <-ggplot(bic_zinph, aes(min_model)) +
  geom_bar(aes(y = (after_stat(count))/sum(after_stat(count)), fill =min_model), show.legend = FALSE) +
  theme_pubclean()+
  labs(x = "Model", y = "% smallest BIC")+ 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))+
  scale_fill_manual(values=c("Poisson" = "#1B9E77", "Quasi" ="#D95F02", 
                             "NB-lin"= "#7570B3", "NB-quad" = "#E7298A", 
                             "ZIP" = "#084594", "ZINB" = "#A6761D" ))+
  coord_fixed(ratio = 2.2)
p7
# finalplot <- ggarrange(p1, p6, p7, ncol = 3, nrow = 1)
# annotate_figure(finalplot, top = text_grob("Parent model:Negbin-quad(strong overdispersion)", 
#                                            color = "red", face = "bold", size = 14))
ggarrange(p1,                                                 # First row with scatter plot
          ggarrange(p6, p7, nrow = 2, labels = c("b", "c")), # Second row with box and dot plots
          ncol = 2, 
          labels = "a"                                        # Labels of the scatter plot
) 


