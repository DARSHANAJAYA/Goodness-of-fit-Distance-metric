# Clear the workspace and set seed for reproducibility
rm(list = ls(all.names = TRUE))
set.seed(123)

# Load necessary libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)
library(ggpubr)

# Set theme for ggplot2
theme_set(theme_pubr())

# Read data files for different models
datafinalpois13 <- read.csv("/users/Shares/djayakumari/poisson/datafinalupdatedpois13.csv")
datafinalpois4648 <- read.csv("/users/Shares/djayakumari/poisson/datafinalupdatedpois4648.csv")
datafinalnegbin13 <- read.csv("/users/Shares/djayakumari/poisson/datafinalupdatednbquad13.csv")
datafinalnegbin4648 <- read.csv("/users/Shares/djayakumari/poisson/datafinalupdatednbquad4648.csv")
datafinalnblin13 <- read.csv("/users/Shares/djayakumari/poisson/datafinalupdatednblin13.csv")
datafinalnblin4648 <- read.csv("/users/Shares/djayakumari/poisson/datafinalupdatednblin4648.csv")
datafinalquasi13 <- read.csv("/users/Shares/djayakumari/poisson/datafinalupdatedquasi13.csv")
datafinalquasi4648 <- read.csv("/users/Shares/djayakumari/poisson/datafinalupdatedquasi4648.csv")

# Merge data for each model
datafinalpois <- rbind(datafinalpois13 %>% slice(rep(1:n(), 15)), 
                       datafinalpois4648 %>% slice(rep(1:n(), 15)))
datafinalnegquad <- rbind(datafinalnegbin13 %>% slice(rep(1:n(), 15)), 
                          datafinalnegbin4648 %>% slice(rep(1:n(), 15)))
datafinalnblin <- rbind(datafinalnblin13 %>% slice(rep(1:n(), 15)), 
                        datafinalnblin4648 %>% slice(rep(1:n(), 15)))
datafinalquasi <- rbind(datafinalquasi13 %>% slice(rep(1:n(), 15)), 
                        datafinalquasi4648 %>% slice(rep(1:n(), 15)))

# Create data frames for each model with required variables
datafinalpois1 <- data.frame(base_value = datafinalpois$base_value_poisson, model1 = "Poisson")
datafinalnegquad1 <- data.frame(base_value = datafinalnegquad$base_value_negbinquad, model1 = "NB-quad")
datafinalnblin1 <- data.frame(base_value = datafinalnblin$base_value_negbinlin, model1 = "NB-lin")
datafinalquasi1 <- data.frame(base_value = datafinalquasi$base_value_quasi, model1 = "Quasi")

# Merge data frames for overall distance plot
overall_dist <- rbind(datafinalpois1, datafinalquasi1, datafinalnblin1, datafinalnegquad1)

# Create a data frame for base distances
base_distance <- data.frame(Poisson = datafinalpois$base_value_poisson,
                            Quasi = datafinalquasi$base_value_quasi,
                            `NB-lin` = datafinalnblin$base_value_negbinlin,
                            `NB-quad` = datafinalnegquad$base_value_negbinquad)

# Add L_norm and sample_size variables to base_distance data frame
L2 <- as.data.frame(rep("L2 norm", times = 44865))
colnames(L2) <- "L"
L1 <- as.data.frame(rep("L1 norm", times = 44865))
colnames(L1) <- "L"
L_norm <- as.data.frame(rbind(L2, L1))

twenty <- rep(20, times = 997)
fifty <- rep(50, times = 997)
hundred <- rep(100, times = 997)
sample_size <- rep(c(twenty, fifty, hundred), times = 30)
sample_size <- as.data.frame(sample_size)

base_distance <- cbind(base_distance, L_norm, sample_size)
base_distance$sample_size <- as.factor(base_distance$sample_size)
base_distance$L <- as.factor(base_distance$L)

# Create a data frame for BIC values
bic_poisson <- data.frame(bic_poisson = datafinalpois$bic_poisson,
                          bic_negbinlin = datafinalnblin$bic_nblin,
                          bic_negbinquad = datafinalnegquad$bic_nbquad,
                          L_norm,
                          sample_size)

bic_poisson$sample_size <- as.factor(bic_poisson$sample_size)
bic_poisson$L <- as.factor(bic_poisson$L)

# Plotting

# Plot for overall distance
p1 <- overall_dist %>% 
  ggplot(aes(x = sample_size, y = log(base_value))) +
  theme_bw() +
  geom_boxplot(aes(color = model1)) +
  facet_wrap(~L) + 
  labs(x = "Sample size", y = "ln(d)") +
  ggtitle(expression("Parent model: Poisson")) +
  theme(legend.position = "top") +
  scale_color_manual(name = "Fitted models",
                     values = c("Poisson" = "#1B9E77", "Quasi" = "#D95F02",
                                "NB-lin" = "#7570B3", "NB-quad" = "#E7298A"))

# Plot for frequency of distances
p6 <- ggplot(base_distance, aes(min_model)) +
  geom_bar(aes(y = (..count..)/sum(..count..), fill = min_model), show.legend = FALSE) +
  theme_pubclean() +
  facet_wrap(~sample_size + L, nrow = 4, labeller = labeller(
    sample_size = ~ paste("sample size:", .),
    L = ~ paste(.),
    .multi_line = FALSE
  )) +
  labs(x = "Model", y = "% smallest d") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = c("Poisson" = "#1B9E77", "Quasi" = "#D95F02",
                               "NB-lin" = "#7570B3", "NB-quad" = "#E7298A"))

# Plot for BIC values
p7 <- ggplot(bic_poisson, aes(min_model)) +
  geom_bar(aes(y = (..count..)/sum(..count..), fill = min_model), show.legend = FALSE) +
  theme_pubclean() +
  facet_wrap(~sample_size + L, nrow = 4, labeller = labeller(
    sample_size = ~ paste("sample size: ", .),
    L = ~ paste(.),
    .multi_line = FALSE
  )) +
  labs(x = "Model", y = "% smallest BIC") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = c("Poisson" = "#1B9E77", "Quasi" = "#D95F02",
                               "NB-lin" = "#7570B3", "NB-quad" = "#E7298A"))

# Arrange plots
ggarrange(p1, ggarrange(p6, p7, nrow = 2, labels = c("b", "c")), ncol = 2, labels = "a")


