#####negbinh

###loading the files 
rm(list = ls(all.names = TRUE))
set.seed(223)

library(plyr)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(reshape2)

library(wesanderson)
library(ggtext)
###initialising the variables 
base_value_poisson <- matrix(nrow = 90, ncol = 498)
base_value_negbinquad <- matrix(nrow = 90, ncol =498)
base_value_negbinlin  <- matrix(nrow = 90, ncol = 498)
base_value_quasi   <- matrix(nrow = 90, ncol =  498)
base_value_zinbh   <- matrix(nrow = 90, ncol = 498)
base_value_zinp <- matrix(nrow = 90, ncol = 498)

b_value_poisson <- matrix(nrow = 90, ncol = 498)
b_value_negbinquad <- matrix(nrow = 90, ncol = 498)
b_value_negbinlin  <- matrix(nrow = 90, ncol = 498)
b_value_quasi   <- matrix(nrow = 90, ncol = 498)
b_value_zinbh <- matrix(nrow = 90, ncol =498)
b_value_zinp <- matrix(nrow = 90, ncol =498)
w_value_poisson <- matrix(nrow = 90, ncol = 498)
w_value_negbinquad <- matrix(nrow = 90, ncol = 498)
w_value_negbinlin  <- matrix(nrow = 90, ncol = 498)
w_value_quasi   <- matrix(nrow = 90, ncol = 498)
w_value_zinbh   <- matrix(nrow = 90, ncol = 498)
w_value_zinp   <- matrix(nrow = 90, ncol = 498)
dist_nb_lin <- matrix(nrow = 90, ncol = 498)
dist_nb_quad <- matrix(nrow = 90, ncol =498)
dist_poisson <- matrix(nrow = 90, ncol = 498)
dist_quasi <- matrix(nrow = 90, ncol = 498)
dist_zinbh <- matrix(nrow = 90, ncol = 498)
dist_zinp <- matrix(nrow = 90, ncol = 498)

dist_w_nb_lin <- matrix(nrow = 90, ncol = 498)
dist_w_nb_quad <- matrix(nrow = 90, ncol = 498)
dist_w_poisson <- matrix(nrow = 90, ncol = 498)
dist_w_quasi <- matrix(nrow = 90, ncol = 498)
dist_w_zinbh <- matrix(nrow = 90, ncol = 498)
dist_w_zinp <- matrix(nrow = 90, ncol = 498)


b_poisson <- matrix(nrow = 90, ncol = 498)
b_nb_quad <- matrix(nrow = 90, ncol = 498)
b_nb_lin  <- matrix(nrow = 90, ncol = 498)
b_quasi   <- matrix(nrow = 90, ncol = 498)
b_zinbh  <- matrix(nrow = 90, ncol = 498)
b_zinp  <- matrix(nrow = 90, ncol = 498)
w_poisson <- matrix(nrow = 90, ncol = 498)
w_nb_quad <- matrix(nrow = 90, ncol = 498)
w_nb_lin  <- matrix(nrow = 90, ncol = 498)
w_quasi   <- matrix(nrow = 90, ncol = 498)
w_zinbh   <- matrix(nrow = 90, ncol = 498)
w_zinp   <- matrix(nrow = 90, ncol = 498)


r_poisson <- matrix(nrow = 90, ncol = 498)
r_nb_quad <- matrix(nrow = 90, ncol = 498)
r_nb_lin  <- matrix(nrow = 90, ncol = 498)
r_quasi   <- matrix(nrow = 90, ncol = 498)
r_zinbh   <- matrix(nrow = 90, ncol = 498)
r_zinp   <- matrix(nrow = 90, ncol = 498)


m_poisson <- matrix(nrow = 90, ncol = 498)
m_nb_quad <- matrix(nrow = 90, ncol = 498)
m_nb_lin  <- matrix(nrow = 90, ncol = 498)
m_quasi   <- matrix(nrow = 90, ncol = 498)
m_zinbh   <- matrix(nrow = 90, ncol = 498)
m_zinp   <- matrix(nrow = 90, ncol = 498)


bic_poisson <- matrix(nrow = 90, ncol = 498)
bic_neg_quad <- matrix(nrow = 90, ncol = 498)
bic_neg_lin <- matrix(nrow =90, ncol = 498)
bic_zinbh <- matrix(nrow =90, ncol = 498)
bic_zinp <- matrix(nrow =90, ncol = 498)

get_poisson       <- list()


#function to load the files 

get_file <- function(file){
  
  load(file)
  data <- my_simulation1
  rm(my_simulation1)
  return(data)
  
  
  
}

#function to manipulate the Rdata files to get the required 

for(j in 1:3) {
  
  
  #setwd("")
  
  get_file_poisson <- paste0("sim_zinbom46_", j , "_new.Rdata")
  
  get_poisson[[j]]<- get_file(get_file_poisson)
  
  for(i in 1:998){  
    # browser()
    # get_poisson[[j]][["w_poisson"]][[i]] <- ifelse(is.null(get_poisson[[j]][["w_poisson"]][[i]]) == TRUE || is.double(get_poisson[[j]][["w_poisson"]][[i]]) == FALSE, 0,get_poisson[[j]][["w_poisson"]][[i]])
    # get_poisson[[j]][["w_nb_quad"]][[i]] <- ifelse(is.null(get_poisson[[j]][["w_nb_quad"]][[i]]) == TRUE || is.double(get_poisson[[j]][["w_nb_quad"]][[i]]) == FALSE, 0,get_poisson[[j]][["w_nb_quad"]][[i]])
    # get_poisson[[j]][["w_nb_lin"]][[i]] <- ifelse(is.null(get_poisson[[j]][["w_nb_lin"]][[i]]) == TRUE || is.double(get_poisson[[j]][["w_nb_lin"]][[i]]) == FALSE, 0,get_poisson[[j]][["w_nb_lin"]][[i]])
    # get_poisson[[j]][["w_quasi"]][[i]] <- ifelse(is.null(get_poisson[[j]][["w_quasi"]][[i]]) == TRUE || is.double(get_poisson[[j]][["w_quasi"]][[i]]) == FALSE, 0,get_poisson[[j]][["w_quasi"]][[i]])
    # #   # d <- replace(df$x, df$x > 4, 5)
    # dist_poisson[j,i]<- get_poisson[[j]][["dist_poisson"]][[i]]
    # # base_dist_poisson[j,i]<- get_poisson[[j]][["dist_poisson"]][[i]]*get_poisson[[j]][["w_poisson"]][[i]]
    # dist_nb_quad[j,i] <- get_poisson[[j]][["dist_nb_quad"]][[i]]
    # 
    # dist_nb_lin[j,i]  <- get_poisson[[j]][["dist_nb_lin"]][[i]]
    # dist_quasi[j,i]  <- get_poisson[[j]][["dist_quasi"]][[i]]
    # dist_zinbh[j,i]  <- get_poisson[[j]][["dist_zinbh"]][[i]]
    # dist_zinp[j,i]  <- get_poisson[[j]][["dist_zinp"]][[i]]

    
    dist_w_poisson[j,i]<- get_poisson[[j]][["dist_w_poisson"]][[i]]
    # base_dist_poisson[j,i]<- get_poisson[[j]][["dist_poisson"]][[i]]*get_poisson[[j]][["w_poisson"]][[i]]
    dist_w_nb_quad[j,i] <- get_poisson[[j]][["dist_w_nb_quad"]][[i]]
    
    dist_w_nb_lin[j,i]  <- get_poisson[[j]][["dist_w_nb_lin"]][[i]]
    dist_w_quasi[j,i]  <- get_poisson[[j]][["dist_w_quasi"]][[i]]
    dist_w_zinbh[j,i]  <- get_poisson[[j]][["dist_w_zinbh"]][[i]]
    dist_w_zinp[j,i]  <- get_poisson[[j]][["dist_w_zinp"]][[i]]
    get_poisson[[j]][["base_value_poisson"]][[i]] <- ifelse(is.double(get_poisson[[j]][["base_value_poisson"]][[i]]) == FALSE, NA,get_poisson[[j]][["base_value_poisson"]][[i]])
    get_poisson[[j]][["base_value_negbinlin"]][[i]] <- ifelse(is.double(get_poisson[[j]][["base_value_negbinlin"]][[i]]) == FALSE, NA,get_poisson[[j]][["base_value_negbinlin"]][[i]])
    get_poisson[[j]][["base_value_negbinquad"]][[i]] <- ifelse(is.double(get_poisson[[j]][["base_value_negbinquad"]][[i]]) == FALSE, NA,get_poisson[[j]][["base_value_negbinquad"]][[i]])
    get_poisson[[j]][[" base_value_quasi"]][[i]] <- ifelse(is.double(get_poisson[[j]][[" base_value_quasi"]][[i]]) == FALSE, NA,get_poisson[[j]][[" base_value_quasi"]][[i]])
    get_poisson[[j]][[" base_value_zinbh"]][[i]] <- ifelse(is.double(get_poisson[[j]][[" base_value_zinbh"]][[i]]) == FALSE, NA,get_poisson[[j]][[" base_value_zinbh"]][[i]])
    get_poisson[[j]][["base_value_zinp"]][[i]] <- ifelse(is.double(get_poisson[[j]][["base_value_zinp"]][[i]]) == FALSE, NA,get_poisson[[j]][["base_value_zinp"]][[i]])
    
    
    
    
    get_poisson[[j]][["b_value_poisson"]][[i]] <- ifelse(is.double(get_poisson[[j]][["b_value_poisson"]][[i]]) == FALSE, NA,get_poisson[[j]][["b_value_poisson"]][[i]])
    get_poisson[[j]][["b_value_negbinlin "]][[i]] <- ifelse(is.double(get_poisson[[j]][["b_value_negbinlin "]][[i]]) == FALSE, NA,get_poisson[[j]][["b_value_negbinlin "]][[i]])
    get_poisson[[j]][["b_value_negbinquad"]][[i]] <- ifelse(is.double(get_poisson[[j]][["b_value_negbinquad"]][[i]]) == FALSE, NA,get_poisson[[j]][["b_value_negbinquad"]][[i]])
    get_poisson[[j]][["b_value_quasi"]][[i]] <- ifelse(is.double(get_poisson[[j]][["b_value_quasi"]][[i]]) == FALSE, NA,get_poisson[[j]][["b_value_quasi"]][[i]])
    get_poisson[[j]][["b_value_zinbh"]][[i]] <- ifelse(is.double(get_poisson[[j]][["b_value_zinbh"]][[i]]) == FALSE, NA,get_poisson[[j]][["b_value_zinbh"]][[i]])
    get_poisson[[j]][["b_value_zinp"]][[i]] <- ifelse(is.double(get_poisson[[j]][["b_value_zinp"]][[i]]) == FALSE, NA,get_poisson[[j]][["b_value_zinp"]][[i]])

    
    
    
    
    
    get_poisson[[j]][["w_value_poisson"]][[i]] <- ifelse(is.double(get_poisson[[j]][["w_value_poisson"]][[i]]) == FALSE, NA,get_poisson[[j]][["w_value_poisson"]][[i]])
    get_poisson[[j]][["w_value_negbinlin"]][[i]] <- ifelse(is.double(get_poisson[[j]][["w_value_negbinlin"]][[i]]) == FALSE, NA,get_poisson[[j]][["w_value_negbinlin"]][[i]])
    get_poisson[[j]][["w_value_negbinquad"]][[i]] <- ifelse(is.double(get_poisson[[j]][["w_value_negbinquad"]][[i]]) == FALSE, NA,get_poisson[[j]][["w_value_negbinquad"]][[i]])
    get_poisson[[j]][["w_value_quasi"]][[i]] <- ifelse(is.double(get_poisson[[j]][["w_value_quasi"]][[i]]) == FALSE, NA,get_poisson[[j]][["w_value_quasi"]][[i]])
    get_poisson[[j]][["w_value_zinbh"]][[i]] <- ifelse(is.double(get_poisson[[j]][["w_value_zinbh"]][[i]]) == FALSE, NA,get_poisson[[j]][["w_value_zinbh"]][[i]])
    get_poisson[[j]][["w_value_zinp"]][[i]] <- ifelse(is.double(get_poisson[[j]][["w_value_zinp"]][[i]]) == FALSE, NA,get_poisson[[j]][["w_value_zinp"]][[i]])
    
    
    get_poisson[[j]][["b_poisson"]][[i]] <- ifelse(is.double(get_poisson[[j]][["b_poisson"]][[i]]) == FALSE, NA,get_poisson[[j]][["b_poisson"]][[i]])
    get_poisson[[j]][["b_nb_lin"]][[i]] <- ifelse(is.double(get_poisson[[j]][["b_nb_lin"]][[i]]) == FALSE, NA,get_poisson[[j]][["b_nb_lin"]][[i]])
    get_poisson[[j]][["b_nb_quad"]][[i]] <- ifelse(is.double(get_poisson[[j]][["b_nb_quad"]][[i]]) == FALSE, NA,get_poisson[[j]][["b_nb_quad"]][[i]])
    get_poisson[[j]][["b_quasi"]][[i]] <- ifelse(is.double(get_poisson[[j]][["b_quasi"]][[i]]) == FALSE, NA,get_poisson[[j]][["b_quasi"]][[i]])
    get_poisson[[j]][["b_zinbh"]][[i]] <- ifelse(is.double(get_poisson[[j]][["b_zinbh"]][[i]]) == FALSE, NA,get_poisson[[j]][["b_zinbh"]][[i]])
    get_poisson[[j]][["b_zinp"]][[i]] <- ifelse(is.double(get_poisson[[j]][["b_zinp"]][[i]]) == FALSE, NA,get_poisson[[j]][["b_zinp"]][[i]])
    
    
    get_poisson[[j]][["w_poisson"]][[i]] <- ifelse(is.double(get_poisson[[j]][["w_poisson"]][[i]]) == FALSE, NA,get_poisson[[j]][["w_poisson"]][[i]])
    get_poisson[[j]][["w_nb_lin"]][[i]] <- ifelse(is.double(get_poisson[[j]][["w_nb_lin"]][[i]]) == FALSE,NA ,get_poisson[[j]][["w_nb_lin"]][[i]])
    get_poisson[[j]][["w_nb_quad"]][[i]] <- ifelse(is.double(get_poisson[[j]][["w_nb_quad"]][[i]]) == FALSE, NA,get_poisson[[j]][["w_nb_quad"]][[i]])
    get_poisson[[j]][["w_quasi"]][[i]] <- ifelse(is.double(get_poisson[[j]][["w_quasi"]][[i]]) == FALSE, NA,get_poisson[[j]][["w_quasi"]][[i]])
    get_poisson[[j]][["w_zinbh"]][[i]] <- ifelse(is.double(get_poisson[[j]][["w_zinbh"]][[i]]) == FALSE, NA,get_poisson[[j]][["w_zinbh"]][[i]])
    get_poisson[[j]][["w_zinp"]][[i]] <- ifelse(is.double(get_poisson[[j]][["w_zinp"]][[i]]) == FALSE, NA,get_poisson[[j]][["w_zinp"]][[i]])
    
    get_poisson[[j]][["r_poisson"]][[i]] <- ifelse(is.double(get_poisson[[j]][["r_poisson"]][[i]]) == FALSE, NA,get_poisson[[j]][["r_poisson"]][[i]])
    get_poisson[[j]][["r_negbinlin"]][[i]] <- ifelse(is.double(get_poisson[[j]][["r_negbinlin"]][[i]]) == FALSE, NA,get_poisson[[j]][["r_negbinlin"]][[i]])
    get_poisson[[j]][["r_negbinquad"]][[i]] <- ifelse(is.double(get_poisson[[j]][["r_negbinquad"]][[i]]) == FALSE, NA,get_poisson[[j]][["r_negbinquad"]][[i]])
    get_poisson[[j]][["r_quasi"]][[i]] <- ifelse(is.double(get_poisson[[j]][["r_quasi"]][[i]]) == FALSE, NA,get_poisson[[j]][["r_quasi"]][[i]])
    get_poisson[[j]][["r_zinbh"]][[i]] <- ifelse(is.double(get_poisson[[j]][["r_zinbh"]][[i]]) == FALSE, NA,get_poisson[[j]][["r_zinbh"]][[i]])
    get_poisson[[j]][["r_zinp"]][[i]] <- ifelse(is.double(get_poisson[[j]][["r_zinp"]][[i]]) == FALSE, NA,get_poisson[[j]][["r_zinp"]][[i]])
    
    get_poisson[[j]][["m_poisson"]][[i]] <- ifelse(is.double(get_poisson[[j]][["m_poisson"]][[i]]) == FALSE, NA,get_poisson[[j]][["m_poisson"]][[i]])
    get_poisson[[j]][["m_negbinlin"]][[i]] <- ifelse(is.double(get_poisson[[j]][["m_negbinlin"]][[i]]) == FALSE, NA,get_poisson[[j]][["m_negbinlin"]][[i]])
    get_poisson[[j]][["m_negbinquad"]][[i]] <- ifelse(is.double(get_poisson[[j]][["m_negbinquad"]][[i]]) == FALSE, NA,get_poisson[[j]][["m_negbinquad"]][[i]])
    get_poisson[[j]][["m_quasi"]][[i]] <- ifelse(is.double(get_poisson[[j]][["m_quasi"]][[i]]) == FALSE, NA,get_poisson[[j]][["m_quasi"]][[i]])
    get_poisson[[j]][["m_zinbh"]][[i]] <- ifelse(is.double(get_poisson[[j]][["m_zinbh"]][[i]]) == FALSE, NA,get_poisson[[j]][["m_zinbh"]][[i]])
    get_poisson[[j]][["m_zinp"]][[i]] <- ifelse(is.double(get_poisson[[j]][["m_zinp"]][[i]]) == FALSE, NA,get_poisson[[j]][["m_zinp"]][[i]])
   
    
    
    get_poisson[[j]][["bic_poisson"]][[i]] <- ifelse(is.double(get_poisson[[j]][["bic_poisson"]][[i]]) == FALSE, NA,get_poisson[[j]][["bic_poisson"]][[i]])
    get_poisson[[j]][["bic_neg_lin"]][[i]] <- ifelse(is.double(get_poisson[[j]][["bic_neg_lin"]][[i]]) == FALSE, NA,get_poisson[[j]][["bic_neg_lin"]][[i]])
    get_poisson[[j]][["bic_neg_quad"]][[i]] <- ifelse(is.double(get_poisson[[j]][["bic_neg_quad"]][[i]]) == FALSE, NA,get_poisson[[j]][["bic_neg_quad"]][[i]])
    
    get_poisson[[j]][["bic_zinbh"]][[i]] <- ifelse(is.double(get_poisson[[j]][["bic_zinbh"]][[i]]) == FALSE, NA,get_poisson[[j]][["bic_zinbh"]][[i]])
    get_poisson[[j]][["bic_zinp"]][[i]] <- ifelse(is.double(get_poisson[[j]][["bic_zinp"]][[i]]) == FALSE, NA,get_poisson[[j]][["bic_zinp"]][[i]])
    
     # 
    base_value_poisson[[j,i]] <- get_poisson[[j]][["base_value_poisson"]][[i]]
    base_value_negbinlin[[j,i]] <- get_poisson[[j]][["base_value_negbinlin"]][[i]]
    base_value_negbinquad[[j,i]] <- get_poisson[[j]][["base_value_negbinquad"]][[i]]
    base_value_quasi[[j,i]]  <- get_poisson[[j]][[" base_value_quasi"]][[i]]
    base_value_zinbh[[j,i]]  <- get_poisson[[j]][[" base_value_zinbh"]][[i]]
    base_value_zinp[[j,i]]  <- get_poisson[[j]][["base_value_zinp"]][[i]]
    
    
    b_value_poisson[[j,i]] <-  get_poisson[[j]][["b_value_poisson"]][[i]]
    b_value_negbinlin[[j,i]] <-  get_poisson[[j]][["b_value_negbinlin "]][[i]]
    b_value_negbinquad[[j,i]] <-  get_poisson[[j]][["b_value_negbinquad"]][[i]]
    b_value_quasi[[j,i]] <-  get_poisson[[j]][["b_value_quasi"]][[i]]
    b_value_zinbh[[j,i]] <-  get_poisson[[j]][["b_value_zinbh"]][[i]]
    b_value_zinp[[j,i]] <-  get_poisson[[j]][["b_value_zinp"]][[i]]
    
    w_value_poisson[[j,i]] <-  get_poisson[[j]][["w_value_poisson"]][[i]]
    w_value_negbinlin[[j,i]] <-  get_poisson[[j]][["w_value_negbinlin"]][[i]]
    w_value_negbinquad[[j,i]] <-  get_poisson[[j]][["w_value_negbinquad"]][[i]]
    w_value_quasi[[j,i]] <-  get_poisson[[j]][["w_value_quasi"]][[i]]
    w_value_zinbh[[j,i]] <-  get_poisson[[j]][["w_value_zinbh"]][[i]]
    w_value_zinp[[j,i]] <-  get_poisson[[j]][["w_value_zinp"]][[i]]

    b_poisson[[j,i]] <-  get_poisson[[j]][["b_poisson"]][[i]]
    b_nb_lin[[j,i]] <-  get_poisson[[j]][["b_nb_lin"]][[i]]
    b_nb_quad[[j,i]] <-  get_poisson[[j]][["b_nb_quad"]][[i]]
    b_quasi[[j,i]] <-  get_poisson[[j]][["b_quasi"]][[i]]
    b_zinbh[[j,i]] <-  get_poisson[[j]][["b_zinbh"]][[i]]
    b_zinp[[j,i]] <-  get_poisson[[j]][["b_zinp"]][[i]]
    
    
    w_poisson[[j,i]] <-  get_poisson[[j]][["w_poisson"]][[i]]
    w_nb_lin[[j,i]] <-  get_poisson[[j]][["w_nb_lin"]][[i]]
    w_nb_quad[[j,i]] <-  get_poisson[[j]][["w_nb_quad"]][[i]]
    w_quasi[[j,i]] <-  get_poisson[[j]][["w_quasi"]][[i]]
    w_zinbh[[j,i]] <-  get_poisson[[j]][["w_zinbh"]][[i]]
    w_zinp[[j,i]] <-  get_poisson[[j]][["w_zinp"]][[i]]
    
    r_poisson[[j,i]] <-  get_poisson[[j]][["r_poisson"]][[i]]
    r_nb_lin[[j,i]] <-  get_poisson[[j]][["r_negbinlin"]][[i]]
    r_nb_quad[[j,i]] <-  get_poisson[[j]][["r_negbinquad"]][[i]]
    r_quasi[[j,i]] <-  get_poisson[[j]][["r_quasi"]][[i]]
    r_zinbh[[j,i]] <-  get_poisson[[j]][["r_zinbh"]][[i]]
    r_zinp[[j,i]] <-  get_poisson[[j]][["r_zinp"]][[i]]
    
    
    m_poisson[[j,i]] <-  get_poisson[[j]][["m_poisson"]][[i]]
    m_nb_lin[[j,i]] <-  get_poisson[[j]][["m_negbinlin"]][[i]]
    m_nb_quad[[j,i]] <-  get_poisson[[j]][["m_negbinquad"]][[i]]
    m_quasi[[j,i]] <-  get_poisson[[j]][["m_quasi"]][[i]]
    m_zinbh[[j,i]] <-  get_poisson[[j]][["m_zinbh"]][[i]]
    m_zinp[[j,i]] <-  get_poisson[[j]][["m_zinp"]][[i]]
    
    
    bic_poisson[[j,i]] <-  get_poisson[[j]][["bic_poisson"]][[i]]
    bic_neg_lin[[j,i]] <-  get_poisson[[j]][["bic_neg_lin"]][[i]]
    bic_neg_quad[[j,i]] <-  get_poisson[[j]][["bic_neg_quad"]][[i]]
    
    bic_zinbh[[j,i]] <-  get_poisson[[j]][["bic_zinbh"]][[i]]
    bic_zinp[[j,i]] <-  get_poisson[[j]][["bic_zinp"]][[i]]

    # 
    # 
  }
  cat("Run number:", j, "\n")  
}




# getting all the variables
# 
# 
base_value_poisson <-unlist(data.frame(t(base_value_poisson[1:3,])))
base_value_negbinlin  <- unlist(data.frame(t(base_value_negbinlin[1:3,])))
base_value_negbinquad <- unlist(data.frame(t(base_value_negbinquad[1:3,])))
base_value_quasi   <- unlist(data.frame(t(base_value_quasi[1:3,])))
base_value_zinbh   <- unlist(data.frame(t(base_value_zinbh[1:3,])))
base_value_zinp   <- unlist(data.frame(t(base_value_zinp[1:3,])))
b_value_poisson <-unlist(data.frame(t(b_value_poisson[1:3,])))
b_value_negbinquad <- unlist(data.frame(t(b_value_negbinquad[1:3,])))
b_value_negbinlin  <- unlist(data.frame(t(b_value_negbinlin[1:3,])))
b_value_quasi   <- unlist(data.frame(t(b_value_quasi[1:3,])))
b_value_zinbh   <- unlist(data.frame(t(b_value_zinbh[1:3,])))
b_value_zinp   <- unlist(data.frame(t(b_value_zinp[1:3,])))
w_value_poisson <- unlist(data.frame(t(w_value_poisson[1:3,])))
w_value_negbinquad <- unlist(data.frame(t(w_value_negbinquad[1:3,])))
w_value_negbinlin  <- unlist(data.frame(t(w_value_negbinlin[1:3,])))
w_value_quasi   <- unlist(data.frame(t(w_value_quasi[1:3,])))
w_value_zinbh  <- unlist(data.frame(t(w_value_zinbh[1:3,])))
w_value_zinp  <- unlist(data.frame(t(w_value_zinp[1:3,])))


b_poisson <- unlist(data.frame(t(b_poisson[1:3,])))
b_nb_quad <- unlist(data.frame(t(b_nb_quad[1:3,])))
b_nb_lin  <- unlist(data.frame(t(b_nb_lin[1:3,])))
b_quasi   <- unlist(data.frame(t(b_quasi[1:3,])))
b_zinbh   <- unlist(data.frame(t(b_zinbh[1:3,])))
b_zinp   <- unlist(data.frame(t(b_zinp[1:3,])))
w_poisson <- unlist(data.frame(t(w_poisson[1:3,])))
w_nb_quad <- unlist(data.frame(t(w_nb_quad[1:3,])))
w_nb_lin  <- unlist(data.frame(t(w_nb_lin[1:3,])))
w_quasi   <- unlist(data.frame(t(w_quasi[1:3,])))
w_zinbh  <- unlist(data.frame(t(w_zinbh[1:3,])))
w_zinp  <- unlist(data.frame(t(w_zinp[1:3,])))

r_poisson <- unlist(data.frame(t(r_poisson[1:3,])))
r_nb_quad <- unlist(data.frame(t(r_nb_quad[1:3,])))
r_nb_lin  <- unlist(data.frame(t(r_nb_lin[1:3,])))
r_quasi   <- unlist(data.frame(t(r_quasi[1:3,])))
r_zinbh   <- unlist(data.frame(t(r_zinbh[1:3,])))
r_zinp   <- unlist(data.frame(t(r_zinp[1:3,])))

m_poisson <- unlist(data.frame(t(m_poisson[1:3,])))
m_nb_quad <- unlist(data.frame(t(m_nb_quad[1:3,])))
m_nb_lin  <- unlist(data.frame(t(m_nb_lin[1:3,])))
m_quasi   <- unlist(data.frame(t(m_quasi[1:3,])))
m_zinbh   <- unlist(data.frame(t(m_zinbh[1:3,])))
m_zinp   <- unlist(data.frame(t(m_zinp[1:3,])))

bic_poisson <- unlist(data.frame(t(bic_poisson[1:3,])))
bic_neg_quad <- unlist(data.frame(t(bic_neg_quad[1:3,])))
bic_neg_lin  <- unlist(data.frame(t(bic_neg_lin[1:3,])))

bic_zinbh   <- unlist(data.frame(t(bic_zinbh[1:3,])))
bic_zinp   <- unlist(data.frame(t(bic_zinp[1:3,])))

# 
dist_w_nb_lin <- unlist(data.frame(t(dist_w_nb_lin[1:3,])))
dist_w_nb_quad <- unlist(data.frame(t(dist_w_nb_quad[1:3,])))
dist_w_poisson <- unlist(data.frame(t(dist_w_poisson[1:3,])))
dist_w_quasi <- unlist(data.frame(t(dist_w_quasi[1:3,])))
dist_w_zinbh <- unlist(data.frame(t(dist_w_zinbh[1:3,])))
dist_w_zinp <- unlist(data.frame(t(dist_w_zinp[1:3,])))



#datacheck <- data.frame(base_value_poisson, b_value_poisson)


datapoisson<- data.frame(base_value_poisson,b_value_poisson,w_value_poisson, b_poisson, w_poisson, r_poisson, m_poisson, dist_w_poisson, bic_poisson)
datanegbinquad <-data.frame(base_value_negbinquad,b_value_negbinquad,w_value_negbinquad, b_nb_quad, w_nb_quad , r_nb_quad, m_nb_quad, dist_w_nb_quad, bic_neg_quad)
datanegbinlin<-data.frame(base_value_negbinlin,b_value_negbinlin,w_value_negbinlin,b_nb_lin, w_nb_lin, r_nb_lin, m_nb_lin, dist_w_nb_lin, bic_neg_lin)
dataquasi <-data.frame(base_value_quasi,b_value_quasi,w_value_quasi,b_quasi, w_quasi, r_quasi, m_quasi, dist_w_quasi)
datazinbh <-data.frame(base_value_zinbh,b_value_zinbh,w_value_zinbh,b_zinbh, w_zinbh, r_zinbh, m_zinbh, dist_w_zinbh, bic_zinbh)
datazinp <-data.frame(base_value_zinp,b_value_zinp,w_value_zinp,b_zinp, w_zinp, r_zinp, m_zinp, dist_w_zinp, bic_zinp)
# colnames(datapoisson22) <- c("base_value_poisson","b_value_poisson","w_value_poisson", "dist_poisson","b_poisson", "w_poisson")


write.csv(datapoisson,"datafinalupdatedpois.csv", row.names = FALSE)
write.csv(datanegbinquad,"datafinalupdatednbquad.csv", row.names = FALSE)
write.csv(datanegbinlin,"datafinalupdatednblin.csv", row.names = FALSE)
write.csv(dataquasi,"datafinalupdatedquasi.csv", row.names = FALSE)
write.csv(datazinbh,"datafinalupdatedzeroinflnb.csv", row.names = FALSE)
write.csv(datazinp,"datafinalupdatedzeroinflp.csv", row.names = FALSE)


