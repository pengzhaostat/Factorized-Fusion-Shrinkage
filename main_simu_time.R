
rm(list=ls())
set.seed(2024)

library(MASS)
library(mvtnorm)
library(foreach)
library(doParallel)

library(pROC)
library(RcppEigen)
library(RcppArmadillo)
library(Rcpp)

file_location = dirname(rstudioapi::getActiveDocumentContext()$path)


setwd(paste(file_location,'/function',sep=""))
source('time_simu_function.R')

cores=detectCores()
cl <- makeCluster(cores[1]) 
registerDoParallel(cl)






# Run simulation
for (n in c(50,100,150,200,250,300,400,500)){
  for (prob in c(0.99)){
    
    score = NULL
    
    # score_option1 =foreach (iter =1:50, .combine='c') %dopar% {
    #   simu_cluster(n = n, k0=2, prob =prob)
    # }
    
    cat('1,\n')
    
    # score_option2 =foreach (iter =1:50, .combine='c') %dopar% {
    #   simu_cluster(n = n, k0=5, prob =prob)
    # }
    
    cat('2,\n')
    
    
    score_option3 =foreach (iter =1:50, .combine='c') %dopar% {
      simu_cluster(n = n, k0=10, prob =prob)
    }
    
    cat('3,\n')
    
    
  #  score = c(score,score_option1,score_option2,score_option3)
    
    score = c(score,score_option3)
    
    
    simu_name = paste(paste(paste("score_time",n,sep = ""),prob,sep = ""),'.csv',sep = "")
    
    write.csv(score, file = simu_name)
  }
}

stopCluster(cl)

