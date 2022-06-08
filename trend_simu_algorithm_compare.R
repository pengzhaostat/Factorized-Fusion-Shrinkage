rm(list=ls())
set.seed(2021)


library(genlasso)
library(MASS)
library(ggplot2)
library(gridExtra)
library(grid)

file_location = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(paste(file_location,'/function',sep=""))
source('MP_gibbs_multi_Sigma_adaptive.R')
source('proximal_gradient.R')
setwd(file_location)



trend_generate = function(T,d){
  X = matrix(rep(0,T*d),nrow = T)
  initial_mean <- rnorm(d)
  X[1,] = initial_mean
  for (t in 2:T) {
    increment = rep(sample(c(0,2),prob = c(0.99,0.01),size=1),d)
    X[t,] <- X[t-1,]+increment*sign(runif(d,-1,1))
  }
  return(X)
}

T=500
d=1

 sigma <- 0.3

Y_true = trend_generate(T,d)

  Y = Y_true + rnorm(T,mean=0, sd =sigma)

   Y =as.matrix(Y,nrow =T)

  d = dim(Y)[2]

par(mfrow=c(1,1))
Sigma_list = vector("list", T)
for (t in 1:T){
  Sigma_list[[t]] = sigma^(-2)*diag(d)  
}


ymax = max(Y+0.5)
ymin = min(Y-0.5)



 MP_mult_list_Sigma = MP_gibbs_mult_Sigma(Y,init_sigma_X1=sigma, Sigma_list=Sigma_list,tau =1000, max_iter =1000,
                                          gap=1e-6,X_true = Y_true,global_prior = 'fixed')


tree_penalized =  function (lambda,tau){
  n = length(lambda)
  sigma_target = diag(rep(1,n))
  beta_target = rep(0,n)
  beta_target[1] = rnorm(1,0,lambda[1]*tau)
  sigma_target[1,1] = lambda[1]^(-1) * tau^(-1)
  for (i in 2:n){
    beta_target[i] = rnorm(1,beta_target[i-1],tau*lambda[i])
    sigma_target[i,i] =lambda[i]^(-1)* tau^(-1)
    sigma_target[i,i-1]=  -lambda[i]^(-1) * tau^(-1)
  }
  sigma_target = solve(sigma_target)
  return(list(beta = beta_target, sigma = sigma_target))
}
K = tree_penalized(rep(1,T),1)$sigma


mylasso <- proxGD(X=K,gamma=0.008, Y=Y,lambda=0.01)

mylasso2 = proxACCE(X=K,gamma=0.006, Y=Y,lambda=0.01)


iter = min(ncol(mylasso$beta),ncol(mylasso2$beta))

err_e_gd = rep(0,iter)
err_p_gd = rep(0,iter)
err_e_agd = rep(0,iter)
err_p_agd = rep(0,iter)

for (i in 1:iter){
  err_e_gd[i] = sqrt(sum((K%*%mylasso$beta[,i]-Y_true)^2/T))
  err_p_gd[i] = sqrt(sum((K%*%mylasso$beta[,i]-Y)^2/T))
  err_e_agd[i] = sqrt(sum((K%*%mylasso2$beta[,i]-Y_true)^2/T))
  err_p_agd[i] = sqrt(sum((K%*%mylasso2$beta[,i]-Y)^2/T))
}  




df_plot = data.frame(list(log_iter=rep(log(1:iter),3),log_estimation_error =c(log(err_e_gd[1:iter]),
                                      log(err_e_agd[1:iter]), log(MP_mult_list_Sigma$err_e[1:iter])), 
                          log_prediction_error = c(log(err_p_gd[1:iter]),
                                                   log(err_p_agd[1:iter]), log(MP_mult_list_Sigma$err[1:iter])),
                          method=c(rep('ISTA',iter),rep('FISTA',iter),rep('MP',iter))))
df_plot[df_plot==-Inf]=NA

p1=ggplot(data=df_plot, aes(x=log_iter, y=log_estimation_error, group=method)) + theme(legend.position="none")+
  geom_line(size=1.5,aes(color = method, linetype = method)) +
  ggtitle("Estimation Convergence") + theme(plot.title = element_text(hjust = 0.5))

p2 = ggplot(data=df_plot, aes(x=log_iter, y=log_prediction_error, group=method)) + theme(legend.position="right")+
  geom_line(size=1.5,aes(color = method, linetype = method)) +
  ggtitle("Prediction Convergence") + theme(plot.title = element_text(hjust = 0.5))

grid.arrange(p1, p2, nrow = 1,widths = c(0.42,0.58))


