rm(list = ls(all.names = TRUE))
set.seed(456)

library(plyr)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(reshape2)


get_poissonnormal <- function(file){
  load(file)
  data_poisson <- my_simulation1
  rm(my_simulation1)
  #data_aic <- list(data_poisson[["aic_poisson_normal"]])
  #data_aic<-as.data.frame(data_aic)
  data_distance<- list(data_poisson[["dist_poisson_normal"]])
  data_distance <- as.data.frame(data_distance)
  #my_list <- cbind(data_aic, data_distance)
  return(data_distance)
}

get_file<-function(filename){
  
  load(filename)
  data<-my_simulation1
  rm(my_simulation1)
  ind1 <- which(sapply(data$model_nb_lin, function(x) class(x)[1] == "try-error"))
  ind2 <- which(sapply(data$model_nb_quad, function(x) class(x)[1] == "try-error"))
  ind3 <- which(sapply(data$model_nb_lin, function(x) x["converged"] == FALSE))
  
  ind4 <- sort(as.numeric(c(ind1,ind2,ind3)))
  
  temporary1 <- data$model_nb_lin
  temporary1[ind4] <- NA
  
  final_indices <- which(sapply(temporary1, function(x) is.na(x) || x$sigma.fv[1] > 1000))
  df<- cbind(data[["dist_poisson"]],data[["dist_nb_lin"]],data[["dist_nb_quad"]], data[["dist_quasi"]] )
  my_list <- list(df, final_indices)
  return(my_list)
}


get_other <- list()
get_poisson <-list()
indices_na <- list()
rm_na <- list()
unlisted <- list()
unlisted_na <- list()
df_final <- list()
mincol <- list()
df_final_dist <-list()
table_dist <- list()
table_dist_final <- list()
table_aic<- list()
get_b <- function(x) {
  
  the_sequence <- c(1:3,16:18,31:33,46:48,61:63,76:78)
  if(x %in% the_sequence) {
    return("one")
  }
  if(x %in% (the_sequence + 3)) {
    return("ratio")
  }
  if(x %in% (the_sequence + 6)) {
    return("hyper")
  }
  if(x %in% (the_sequence + 9)) {
    return("logistic")
  }
  if(x %in% (the_sequence + 12)) {
    return("linear")
  }
}
get_w<-function(x) {
  the_wseq<- c(1:15,46:60)
  if(x %in% the_wseq) {
    return('one')
  }
  if(x %in% (the_wseq + 15)) {
    return('w')
  }
  if(x %in% (the_wseq + 30)) {
    return('w^2')
  }
}

get_distancenorm <-function(x) {
  the_lseq <- 1:45
  if(x %in% the_lseq){
    return('L2')
  }
  if(x %in% (the_lseq + 45)){
    return('L1')
  }
}
get_samplesize<-function(x) {
  the_sampleseq<-seq(1,90,3)
  if(x %in% the_sampleseq) {
    return(20)
  }
  if(x %in% (the_sampleseq + 1)){
    return(50)
  }
  if(x %in% (the_sampleseq+2)){
    return(100)
  }
}
b <- character(0)
function_w <- character(0)
distance_norm <- character(0)
sample_size <- numeric(0)


for(i in 1:90) {
  if(i %in% 1:45) {
    setwd('/Volumes/Darshana_PhD/project/l2norm/poissonnormal/negbinhlin')
  } else if(i %in% 46:90){
    setwd('/Users/darshanaj/Documents/individual project/sim_new/L1 norm/poissonnormal/negbinhlin')
    
  }
  get_file_poisson <- paste0("sim_negbinhlin_poisson_normal", i, ".RData")
  
  get_poisson[[i]]<- get_poissonnormal(get_file_poisson)
  
  if(i %in% 1:45) {
    setwd('/Volumes/Darshana_PhD/project/L2 norm_new/negbinhlin')
  } else if(i %in% 46:90){
    setwd('/Volumes/Darshana_PhD/project/L1 norm_new/negbinhlin')
    
  }
  
  filename <- paste0("sim_negbinh_lin", i,".RData")
  
  get_other[[i]]<- get_file(filename)
  indices_na[[i]] <-get_other[[i]][[2]]
  get_poisson[[i]][[1]][indices_na[[i]]] <- NA
  unlisted[[i]] <- unlist(get_poisson[[i]])
  unlisted_na[[i]] <- na.omit(unlisted[[i]])
  df_final[[i]] <- as.data.frame(cbind(get_other[[i]][[1]],unlisted_na[[i]]))
  colnames(df_final[[i]]) <- c("dist_poisson","dist_nb_lin","dist_nb_quad","dist_quasi","dist_poisson_normal")
  mincol[[i]] <- colnames(df_final[[i]])[apply(df_final[[i]],1,which.min)]
  #mincol[[i]] <- data.frame(mincol[[i]]) 
  df_final_dist[[i]] <- cbind(df_final[[i]], mincol[[i]])
  
  table_aic[[i]]<-data.frame(table(df_final_dist[[i]]$mincol))
  table_dist[[i]]<-df_final_dist[[i]] %>% count(freq= factor(df_final_dist[[i]]$mincol, c("dist_nb_lin","dist_nb_quad","dist_poisson","dist_poisson_normal","dist_quasi")), .drop = FALSE)
  table_dist_final[[i]] <- table_dist[[i]][[2]]
  b[i] <- get_b(i)
  function_w[i]<-get_w(i)
  distance_norm[i]<- get_distancenorm(i)
  sample_size[i] <- get_samplesize(i)
  
  
  
  
  
  cat(i,"\n")
}


temp<-ldply(table_dist_final, rbind)
table_dist_final <-data.frame(temp)
colnames(table_dist_final) <- c("dist_nb_lin","dist_nb_quad","dist_poisson","dist_poisson_normal","dist_quasi")
table_dist_final$b <- b
#df_distance$b <- b
table_dist_final$function_w <- function_w
#df_distance$function_w <- function_w
table_dist_final$distance_norm <- distance_norm
#df_distance$distance_norm <- distance_norm
table_dist_final$sample_size <- sample_size
#df_distance$sample_size <- sample_size
write.csv(table_dist_final,"/Users/darshanaj/Documents/individual project/sim_new/dataframes/negbinhlin_mindist.csv", row.names = FALSE)



full_pre_processed <- list()
singledata_aic<-list()
temp<-list()
temp1<-list()
singledata_dist<-list()

get_b <- function(x) {
  
  the_sequence <- c(1:3,16:18,31:33,46:48,61:63,76:78)
  if(x %in% the_sequence) {
    return("one")
  }
  if(x %in% (the_sequence + 3)) {
    return("ratio")
  }
  if(x %in% (the_sequence + 6)) {
    return("hyper")
  }
  if(x %in% (the_sequence + 9)) {
    return("logistic")
  }
  if(x %in% (the_sequence + 12)) {
    return("linear")
  }
}
get_w<-function(x) {
  the_wseq<- c(1:15,46:60)
  if(x %in% the_wseq) {
    return('one')
  }
  if(x %in% (the_wseq + 15)) {
    return('w')
  }
  if(x %in% (the_wseq + 30)) {
    return('w^2')
  }
}

get_distancenorm <-function(x) {
  the_lseq <- 1:45
  if(x %in% the_lseq){
    return('L2')
  }
  if(x %in% (the_lseq + 45)){
    return('L1')
  }
}
get_samplesize<-function(x) {
  the_sampleseq<-seq(1,90,3)
  if(x %in% the_sampleseq) {
    return(20)
  }
  if(x %in% (the_sampleseq + 1)){
    return(50)
  }
  if(x %in% (the_sampleseq+2)){
    return(100)
  }
}
