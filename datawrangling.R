#####negbinh


rm(list = ls(all.names = TRUE))
set.seed(6248)
##loadng the libraries 
library(plyr)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(reshape2)

library(wesanderson)
library(ggtext)
###initialising the variables 
base_value_poisson <- matrix(nrow = 90, ncol = 997)
base_value_negbinquad <- matrix(nrow = 90, ncol = 997)
base_value_negbinlin  <- matrix(nrow = 90, ncol = 997)
base_value_quasi   <- matrix(nrow = 90, ncol = 997)
base_value_zinbh   <- matrix(nrow = 90, ncol = 997)
base_value_zinp <- matrix(nrow = 90, ncol = 997)

b_value_poisson <- matrix(nrow = 90, ncol = 997)
b_value_negbinquad <- matrix(nrow = 90, ncol = 997)
b_value_negbinlin  <- matrix(nrow = 90, ncol = 997)
b_value_quasi   <- matrix(nrow = 90, ncol =997)
b_value_zinbh <- matrix(nrow = 90, ncol =997)
b_value_zinp <- matrix(nrow = 90, ncol =997)
w_value_poisson <- matrix(nrow = 90, ncol = 997)
w_value_negbinquad <- matrix(nrow = 90, ncol = 997)
w_value_negbinlin  <- matrix(nrow = 90, ncol = 997)
w_value_quasi   <- matrix(nrow = 90, ncol = 997)
w_value_zinbh   <- matrix(nrow = 90, ncol = 997)
w_value_zinp   <- matrix(nrow = 90, ncol = 997)
dist_nb_lin <- matrix(nrow = 90, ncol = 997)
dist_nb_quad <- matrix(nrow = 90, ncol = 997)
dist_poisson <- matrix(nrow = 90, ncol = 997)
dist_quasi <- matrix(nrow = 90, ncol = 997)
dist_zinbh <- matrix(nrow = 90, ncol = 997)
dist_zinp <- matrix(nrow = 90, ncol = 997)

dist_w_nb_lin <- matrix(nrow = 90, ncol = 997)
dist_w_nb_quad <- matrix(nrow = 90, ncol = 997)
dist_w_poisson <- matrix(nrow = 90, ncol = 997)
dist_w_quasi <- matrix(nrow = 90, ncol = 997)
dist_w_zinbh <- matrix(nrow = 90, ncol = 997)
dist_w_zinp <- matrix(nrow = 90, ncol = 997)


b_poisson <- matrix(nrow = 90, ncol = 997)
b_nb_quad <- matrix(nrow = 90, ncol = 997)
b_nb_lin  <- matrix(nrow = 90, ncol = 997)
b_quasi   <- matrix(nrow = 90, ncol = 997)
b_zinbh  <- matrix(nrow = 90, ncol = 997)
b_zinp  <- matrix(nrow = 90, ncol =997)
w_poisson <- matrix(nrow = 90, ncol = 997)
w_nb_quad <- matrix(nrow = 90, ncol = 997)
w_nb_lin  <- matrix(nrow = 90, ncol = 997)
w_quasi   <- matrix(nrow = 90, ncol = 997)
w_zinbh   <- matrix(nrow = 90, ncol = 997)
w_zinp   <- matrix(nrow = 90, ncol = 997)


r_poisson <- matrix(nrow = 90, ncol = 997)
r_nb_quad <- matrix(nrow = 90, ncol = 997)
r_nb_lin  <- matrix(nrow = 90, ncol = 997)
r_quasi   <- matrix(nrow = 90, ncol = 997)
r_zinbh   <- matrix(nrow = 90, ncol = 997)
r_zinp   <- matrix(nrow = 90, ncol = 997)


m_poisson <- matrix(nrow = 90, ncol = 997)
m_nb_quad <- matrix(nrow = 90, ncol = 997)
m_nb_lin  <- matrix(nrow = 90, ncol = 997)
m_quasi   <- matrix(nrow = 90, ncol = 997)
m_zinbh   <- matrix(nrow = 90, ncol = 997)
m_zinp   <- matrix(nrow = 90, ncol = 997)


bic_poisson <- matrix(nrow = 90, ncol = 997)
bic_nbquad <- matrix(nrow = 90, ncol = 997)
bic_nblin <- matrix(nrow =90, ncol = 997)

get_poisson       <- list()




###function to load the simulation files 

get_file <- function(file){
  
  load(file)
  data <- my_simulation1
  rm(my_simulation1)
  return(data)
  
  
  
}

###function to extrat the essential variables from the rdata file

for(j in 1:3) {
  
  
  #setwd("")
  
  get_file_poisson <- paste0("sim_pois", j , "_new.Rdata")
  
  get_poisson[[j]]<- get_file(get_file_poisson)
  
  for(i in 1:997){  
   
    
    dist_w_poisson[j,i]<- get_poisson[[j]][["dist_w_poisson"]][[i]]
    
    dist_w_nb_quad[j,i] <- get_poisson[[j]][["dist_w_nb_quad"]][[i]]
    
    dist_w_nb_lin[j,i]  <- get_poisson[[j]][["dist_w_nb_lin"]][[i]]
    dist_w_quasi[j,i]  <- get_poisson[[j]][["dist_w_quasi"]][[i]]
   
    get_poisson[[j]][["base_value_poisson"]][[i]] <- ifelse(is.double(get_poisson[[j]][["base_value_poisson"]][[i]]) == FALSE, NA,get_poisson[[j]][["base_value_poisson"]][[i]])
    get_poisson[[j]][["base_value_negbinlin"]][[i]] <- ifelse(is.double(get_poisson[[j]][["base_value_negbinlin"]][[i]]) == FALSE, NA,get_poisson[[j]][["base_value_negbinlin"]][[i]])
    get_poisson[[j]][["base_value_negbinquad"]][[i]] <- ifelse(is.double(get_poisson[[j]][["base_value_negbinquad"]][[i]]) == FALSE, NA,get_poisson[[j]][["base_value_negbinquad"]][[i]])
    get_poisson[[j]][[" base_value_quasi"]][[i]] <- ifelse(is.double(get_poisson[[j]][[" base_value_quasi"]][[i]]) == FALSE, NA,get_poisson[[j]][[" base_value_quasi"]][[i]])
   
    
    
    
    get_poisson[[j]][["b_value_poisson"]][[i]] <- ifelse(is.double(get_poisson[[j]][["b_value_poisson"]][[i]]) == FALSE, NA,get_poisson[[j]][["b_value_poisson"]][[i]])
    get_poisson[[j]][["b_value_negbinlin "]][[i]] <- ifelse(is.double(get_poisson[[j]][["b_value_negbinlin "]][[i]]) == FALSE, NA,get_poisson[[j]][["b_value_negbinlin "]][[i]])
    get_poisson[[j]][["b_value_negbinquad"]][[i]] <- ifelse(is.double(get_poisson[[j]][["b_value_negbinquad"]][[i]]) == FALSE, NA,get_poisson[[j]][["b_value_negbinquad"]][[i]])
    get_poisson[[j]][["b_value_quasi"]][[i]] <- ifelse(is.double(get_poisson[[j]][["b_value_quasi"]][[i]]) == FALSE, NA,get_poisson[[j]][["b_value_quasi"]][[i]])
   
    
    
    
    
    get_poisson[[j]][["w_value_poisson"]][[i]] <- ifelse(is.double(get_poisson[[j]][["w_value_poisson"]][[i]]) == FALSE, NA,get_poisson[[j]][["w_value_poisson"]][[i]])
    get_poisson[[j]][["w_value_negbinlin"]][[i]] <- ifelse(is.double(get_poisson[[j]][["w_value_negbinlin"]][[i]]) == FALSE, NA,get_poisson[[j]][["w_value_negbinlin"]][[i]])
    get_poisson[[j]][["w_value_negbinquad"]][[i]] <- ifelse(is.double(get_poisson[[j]][["w_value_negbinquad"]][[i]]) == FALSE, NA,get_poisson[[j]][["w_value_negbinquad"]][[i]])
    get_poisson[[j]][["w_value_quasi"]][[i]] <- ifelse(is.double(get_poisson[[j]][["w_value_quasi"]][[i]]) == FALSE, NA,get_poisson[[j]][["w_value_quasi"]][[i]])
    
    get_poisson[[j]][["b_poisson"]][[i]] <- ifelse(is.double(get_poisson[[j]][["b_poisson"]][[i]]) == FALSE, NA,get_poisson[[j]][["b_poisson"]][[i]])
    get_poisson[[j]][["b_nb_lin"]][[i]] <- ifelse(is.double(get_poisson[[j]][["b_nb_lin"]][[i]]) == FALSE, NA,get_poisson[[j]][["b_nb_lin"]][[i]])
    get_poisson[[j]][["b_nb_quad"]][[i]] <- ifelse(is.double(get_poisson[[j]][["b_nb_quad"]][[i]]) == FALSE, NA,get_poisson[[j]][["b_nb_quad"]][[i]])
    get_poisson[[j]][["b_quasi"]][[i]] <- ifelse(is.double(get_poisson[[j]][["b_quasi"]][[i]]) == FALSE, NA,get_poisson[[j]][["b_quasi"]][[i]])
    
    get_poisson[[j]][["w_poisson"]][[i]] <- ifelse(is.double(get_poisson[[j]][["w_poisson"]][[i]]) == FALSE, NA,get_poisson[[j]][["w_poisson"]][[i]])
    get_poisson[[j]][["w_nb_lin"]][[i]] <- ifelse(is.double(get_poisson[[j]][["w_nb_lin"]][[i]]) == FALSE,NA ,get_poisson[[j]][["w_nb_lin"]][[i]])
    get_poisson[[j]][["w_nb_quad"]][[i]] <- ifelse(is.double(get_poisson[[j]][["w_nb_quad"]][[i]]) == FALSE, NA,get_poisson[[j]][["w_nb_quad"]][[i]])
    get_poisson[[j]][["w_quasi"]][[i]] <- ifelse(is.double(get_poisson[[j]][["w_quasi"]][[i]]) == FALSE, NA,get_poisson[[j]][["w_quasi"]][[i]])
    
    get_poisson[[j]][["r_poisson"]][[i]] <- ifelse(is.double(get_poisson[[j]][["r_poisson"]][[i]]) == FALSE, NA,get_poisson[[j]][["r_poisson"]][[i]])
    get_poisson[[j]][["r_negbinlin"]][[i]] <- ifelse(is.double(get_poisson[[j]][["r_negbinlin"]][[i]]) == FALSE, NA,get_poisson[[j]][["r_negbinlin"]][[i]])
    get_poisson[[j]][["r_negbinquad"]][[i]] <- ifelse(is.double(get_poisson[[j]][["r_negbinquad"]][[i]]) == FALSE, NA,get_poisson[[j]][["r_negbinquad"]][[i]])
    get_poisson[[j]][["r_quasi"]][[i]] <- ifelse(is.double(get_poisson[[j]][["r_quasi"]][[i]]) == FALSE, NA,get_poisson[[j]][["r_quasi"]][[i]])
    # get_poisson[[j]][["r_zinbh"]][[i]] <- ifelse(is.double(get_poisson[[j]][["r_zinbh"]][[i]]) == FALSE, NA,get_poisson[[j]][["r_zinbh"]][[i]])
    # get_poisson[[j]][["r_zinp"]][[i]] <- ifelse(is.double(get_poisson[[j]][["r_zinp"]][[i]]) == FALSE, NA,get_poisson[[j]][["r_zinp"]][[i]])
    
    get_poisson[[j]][["m_poisson"]][[i]] <- ifelse(is.double(get_poisson[[j]][["m_poisson"]][[i]]) == FALSE, NA,get_poisson[[j]][["m_poisson"]][[i]])
    get_poisson[[j]][["m_negbinlin"]][[i]] <- ifelse(is.double(get_poisson[[j]][["m_negbinlin"]][[i]]) == FALSE, NA,get_poisson[[j]][["m_negbinlin"]][[i]])
    get_poisson[[j]][["m_negbinquad"]][[i]] <- ifelse(is.double(get_poisson[[j]][["m_negbinquad"]][[i]]) == FALSE, NA,get_poisson[[j]][["m_negbinquad"]][[i]])
    get_poisson[[j]][["m_quasi"]][[i]] <- ifelse(is.double(get_poisson[[j]][["m_quasi"]][[i]]) == FALSE, NA,get_poisson[[j]][["m_quasi"]][[i]])
   
    
    get_poisson[[j]][["bic_poisson"]][[i]] <- ifelse(is.double(get_poisson[[j]][["bic_poisson"]][[i]]) == FALSE, 0,get_poisson[[j]][["bic_poisson"]][[i]])
    get_poisson[[j]][["bic_neg_lin"]][[i]] <- ifelse(is.double(get_poisson[[j]][["bic_neg_lin"]][[i]]) == FALSE, 0,get_poisson[[j]][["bic_neg_lin"]][[i]])
    get_poisson[[j]][["bic_neg_quad"]][[i]] <- ifelse(is.double( get_poisson[[j]][["bic_neg_quad"]][[i]]) == FALSE, 0, get_poisson[[j]][["bic_neg_quad"]][[i]])
    
   
    bic_poisson[[j,i]] <-  get_poisson[[j]][["bic_poisson"]][[i]]
    bic_nblin[[j,i]] <-  get_poisson[[j]][["bic_neg_lin"]][[i]]
    bic_nbquad[[j,i]] <-  get_poisson[[j]][["bic_neg_quad"]][[i]]
    
    # 
    base_value_poisson[[j,i]] <- get_poisson[[j]][["base_value_poisson"]][[i]]
    base_value_negbinlin[[j,i]] <- get_poisson[[j]][["base_value_negbinlin"]][[i]]
    base_value_negbinquad[[j,i]] <- get_poisson[[j]][["base_value_negbinquad"]][[i]]
    base_value_quasi[[j,i]]  <- get_poisson[[j]][[" base_value_quasi"]][[i]]
   
    
    b_value_poisson[[j,i]] <-  get_poisson[[j]][["b_value_poisson"]][[i]]
    b_value_negbinlin[[j,i]] <-  get_poisson[[j]][["b_value_negbinlin "]][[i]]
    b_value_negbinquad[[j,i]] <-  get_poisson[[j]][["b_value_negbinquad"]][[i]]
    b_value_quasi[[j,i]] <-  get_poisson[[j]][["b_value_quasi"]][[i]]
   
    w_value_poisson[[j,i]] <-  get_poisson[[j]][["w_value_poisson"]][[i]]
    w_value_negbinlin[[j,i]] <-  get_poisson[[j]][["w_value_negbinlin"]][[i]]
    w_value_negbinquad[[j,i]] <-  get_poisson[[j]][["w_value_negbinquad"]][[i]]
    w_value_quasi[[j,i]] <-  get_poisson[[j]][["w_value_quasi"]][[i]]
    
    b_poisson[[j,i]] <-  get_poisson[[j]][["b_poisson"]][[i]]
    b_nb_lin[[j,i]] <-  get_poisson[[j]][["b_nb_lin"]][[i]]
    b_nb_quad[[j,i]] <-  get_poisson[[j]][["b_nb_quad"]][[i]]
    b_quasi[[j,i]] <-  get_poisson[[j]][["b_quasi"]][[i]]
    # b_zinbh[[j,i]] <-  get_poisson[[j]][["b_zinbh"]][[i]]
    # b_zinp[[j,i]] <-  get_poisson[[j]][["b_zinp"]][[i]]
    
    
    w_poisson[[j,i]] <-  get_poisson[[j]][["w_poisson"]][[i]]
    w_nb_lin[[j,i]] <-  get_poisson[[j]][["w_nb_lin"]][[i]]
    w_nb_quad[[j,i]] <-  get_poisson[[j]][["w_nb_quad"]][[i]]
    w_quasi[[j,i]] <-  get_poisson[[j]][["w_quasi"]][[i]]
    
    r_poisson[[j,i]] <-  get_poisson[[j]][["r_poisson"]][[i]]
    r_nb_lin[[j,i]] <-  get_poisson[[j]][["r_negbinlin"]][[i]]
    r_nb_quad[[j,i]] <-  get_poisson[[j]][["r_negbinquad"]][[i]]
    r_quasi[[j,i]] <-  get_poisson[[j]][["r_quasi"]][[i]]
    
    
    m_poisson[[j,i]] <-  get_poisson[[j]][["m_poisson"]][[i]]
    m_nb_lin[[j,i]] <-  get_poisson[[j]][["m_negbinlin"]][[i]]
    m_nb_quad[[j,i]] <-  get_poisson[[j]][["m_negbinquad"]][[i]]
    m_quasi[[j,i]] <-  get_poisson[[j]][["m_quasi"]][[i]]
  
  }
  cat("Run number:", j, "\n")  
}
#subsetting for finding base distances 
base_value_poisson <-unlist(data.frame(t(base_value_poisson[1:3,])))
base_value_negbinlin  <- unlist(data.frame(t(base_value_negbinlin[1:3,])))
base_value_negbinquad <- unlist(data.frame(t(base_value_negbinquad[1:3,])))
base_value_quasi   <- unlist(data.frame(t(base_value_quasi[1:3,])))

b_value_poisson <-unlist(data.frame(t(b_value_poisson[1:3,])))
b_value_negbinquad <- unlist(data.frame(t(b_value_negbinquad[1:3,])))
b_value_negbinlin  <- unlist(data.frame(t(b_value_negbinlin[1:3,])))
b_value_quasi   <- unlist(data.frame(t(b_value_quasi[1:3,])))
w_value_poisson <- unlist(data.frame(t(w_value_poisson[1:3,])))
w_value_negbinquad <- unlist(data.frame(t(w_value_negbinquad[1:3,])))
w_value_negbinlin  <- unlist(data.frame(t(w_value_negbinlin[1:3,])))
w_value_quasi   <- unlist(data.frame(t(w_value_quasi[1:3,])))

b_poisson <- unlist(data.frame(t(b_poisson[1:3,])))
b_nb_quad <- unlist(data.frame(t(b_nb_quad[1:3,])))
b_nb_lin  <- unlist(data.frame(t(b_nb_lin[1:3,])))
b_quasi   <- unlist(data.frame(t(b_quasi[1:3,])))

w_poisson <- unlist(data.frame(t(w_poisson[1:3,])))
w_nb_quad <- unlist(data.frame(t(w_nb_quad[1:3,])))
w_nb_lin  <- unlist(data.frame(t(w_nb_lin[1:3,])))
w_quasi   <- unlist(data.frame(t(w_quasi[1:3,])))


r_poisson <- unlist(data.frame(t(r_poisson[1:3,])))
r_nb_quad <- unlist(data.frame(t(r_nb_quad[1:3,])))
r_nb_lin  <- unlist(data.frame(t(r_nb_lin[1:3,])))
r_quasi   <- unlist(data.frame(t(r_quasi[1:3,])))


m_poisson <- unlist(data.frame(t(m_poisson[1:3,])))
m_nb_quad <- unlist(data.frame(t(m_nb_quad[1:3,])))
m_nb_lin  <- unlist(data.frame(t(m_nb_lin[1:3,])))
m_quasi   <- unlist(data.frame(t(m_quasi[1:3,])))

# 
dist_w_nb_lin <- unlist(data.frame(t(dist_w_nb_lin[1:3,])))
dist_w_nb_quad <- unlist(data.frame(t(dist_w_nb_quad[1:3,])))
dist_w_poisson <- unlist(data.frame(t(dist_w_poisson[1:3,])))
dist_w_quasi <- unlist(data.frame(t(dist_w_quasi[1:3,])))

bic_poisson <- unlist(data.frame(t(bic_poisson[1:3,])))
bic_nbquad <- unlist(data.frame(t(bic_nbquad[1:3,])))
bic_nblin  <- unlist(data.frame(t(bic_nblin[1:3,])))




datapoisson<- data.frame(base_value_poisson,b_value_poisson,w_value_poisson, b_poisson, w_poisson, r_poisson, m_poisson, dist_w_poisson, bic_poisson)
datanegbinquad <-data.frame(base_value_negbinquad,b_value_negbinquad,w_value_negbinquad, b_nb_quad, w_nb_quad , r_nb_quad, m_nb_quad, dist_w_nb_quad, bic_nbquad)
datanegbinlin<-data.frame(base_value_negbinlin,b_value_negbinlin,w_value_negbinlin,b_nb_lin, w_nb_lin, r_nb_lin, m_nb_lin, dist_w_nb_lin, bic_nblin)
dataquasi <-data.frame(base_value_quasi,b_value_quasi,w_value_quasi,b_quasi, w_quasi, r_quasi, m_quasi, dist_w_quasi)

write.csv(datapoisson,"datafinalupdatedpois13.csv", row.names = FALSE)
write.csv(datanegbinquad,"/datafinalupdatednbquad13.csv", row.names = FALSE)
write.csv(datanegbinlin,"datafinalupdatednblin13.csv", row.names = FALSE)
write.csv(dataquasi,"datafinalupdatedquasi13.csv", row.names = FALSE)

=
