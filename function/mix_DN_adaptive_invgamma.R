#----------------------------------------------------------------------------
# Binary network via SMF variational inference
# Input Variables:

# trans_sd: transition sd
# Y: network values, list of length T, Y[[t]]: n*n network values at time point t
# mean_beta_prior,sigma_beta_prior: mean and sd for prior of beta
# rho: ar1 coefficient, fixed to 1 in the paper
# sigma_X1: prior sd for the initial state
# gap: gap between errors for convergence
# max_iter: maximal cycles in computation
# alpha: fractional power of the likelihood, fixed to be 0.95 in this paper

#Output variables:
# train_AUC: dynamic of the training AUC
# Mean_X: variational mean of each node, list of length T, Mean_X[[t]]: n*d matrix at time point t
# Sigma_X: variational covariance matrix of each node 
# mean_beta: variational mean of intercept beta
# sigma_beta: variational sd of intercept beta 
# iter: k-1
# AUC: AUC for final variational mean

mix_DN_adaptive_invgamma = function(Y,mean_beta_prior=0, sigma_beta_prior=sqrt(10),  gap = 0.01 ,min_iter=1 ,rho=1,alpha=0.95, trans_sd_init =0.1){
  
  
  
  T = length(Y)
  
  n = length(Y[[1]][1,])
  
  
  #### target parameters
  
  mean_beta = 0;    #mean of beta, 
  
  sigma_beta = 0.5;  #initial sd of beta, 
  
  Mean_X = vector("list", T); #Mean of X, Mean_X[[t]][i,] 
  
  Forward.Mean = vector("list", T)
  
  Backward.Mean = vector("list", T)
  
  Sigma_X = replicate(n=T, expr=list()); # covariance matrix of X, Sigma_x[[t]][[i]] is the covariance matrix Sigma_{it}
  
  Forward.var = replicate(n=T, expr=list())
  
  Backward.var = replicate(n=T, expr=list())
  
  Cross_cov = replicate(n=T, expr=list())
  
  for (t in 1:T){
    Mean_X[[t]] = matrix(rep(0,n*d),nrow = n)
    Forward.Mean[[t]] = matrix(rep(0,n*d),nrow = n)
    Backward.Mean[[t]] = matrix(rep(0,n*d),nrow = n)
    Sigma_X[[t]] = vector("list", n)
    Forward.var[[t]]= vector("list", n) 
    Backward.var[[t]] = vector("list", n) 
    Cross_cov[[t]] = vector("list", n)
    for (i in 1:n){
      Mean_X[[t]][i,] =  rmvnorm(n=1,mean=rep(0,d), sigma=diag(d)) # randomly initialization for mu_{it}
      Sigma_X[[t]][[i]] =  diag(d)          # initialization for Sigma_{it}
      Forward.var[[t]][[i]] = diag(d) 
      Backward.var[[t]][[i]] = diag(d) 
      Cross_cov[[t]][[i]] = 0*diag(2*d)
    }
  }
  
  
  #### tangent parameters
  
  Xi = vector("list", T)

  
  for (t in 1:T) { 
    Xi[[t]] = matrix(rep(1,n*n),nrow = n)
  }
  
  
  axis_mat = matrix(2*runif((n)*d)-1,nrow = n)
  #--------------------------------Algorithm---------------------
  K= 500
  
  err =rep(0,K)
  
  train_auc =rep(0,K)
  
  train_auc[1]=-10
  
  k = 2
  
  ind = 0
  
  trans_sd = trans_sd_init
  
  sigma_X1 = 0.2
  
  while(k<K && ind ==0){
    
    
    # Calculate auxiliary values
    
    
    V_beta_cumulative = 0
    
    M_beta_cumulative = 0
    
    V_X_cumulative =  replicate(n=T, expr=list())
    
    M_X_cumulative = vector("list", T)
    
    
    for (t in 1:T){
      M_X_cumulative[[t]] = matrix(rep(0,n*d),nrow = n)
      V_X_cumulative[[t]] = vector("list", n)
      for (i in 1:n){
        V_X_cumulative[[t]][[i]] =  0*diag(d)        
      }
    }
    
    
    for (t in 1:T){
      for( i in 1:n){
        for (j in 1:n){
          M_i = Mean_X[[t]][i,]
          M_j = Mean_X[[t]][j,]
          V_i = Sigma_X[[t]][[i]]
          V_j = Sigma_X[[t]][[j]]
          c = sqrt(distance_squared_inner_prod(M_i,M_j,V_i,V_j) +mean_beta^2+sigma_beta^2)
          Xi[[t]][i,j] = 1/(2*c)*(exp(c)-1)/(exp(c)+1)/(-2)
          if (j!= i){ 
            V_beta_cumulative= V_beta_cumulative -  2*Xi[[t]][i,j]*alpha
            M_beta_cumulative = M_beta_cumulative + (Y[[t]][i,j]-0.5+2*Xi[[t]][i,j]* t(M_i) %*% M_j)*alpha
          }
        }
      }
    }
    
    
    
    
    mean_beta_new = 0;    
    
    sigma_beta_new = 0  
    
    Mean_X_new = vector("list", T)
    
    Sigma_X_new = replicate(n=T, expr=list())
    
    Forward.Mean.new = vector("list", T)
    
    Backward.Mean.new = vector("list", T)
    
    Forward.var.new = replicate(n=T, expr=list())
    
    Backward.var.new = replicate(n=T, expr=list())
    
    
    for (t in 1:T){
      Mean_X_new[[t]] = matrix(rep(0,n*d),nrow = n)
      Forward.Mean.new[[t]] = matrix(rep(0,n*d),nrow = n)
      Backward.Mean.new[[t]] = matrix(rep(0,n*d),nrow = n)
      Sigma_X_new[[t]] = vector("list", n)
      Forward.var.new[[t]] = vector("list", n)
      Backward.var.new[[t]] = vector("list", n)
      for (i in 1:n){
        Sigma_X_new[[t]][[i]] =  0*diag(d)       
        Forward.var.new[[t]][[i]] =  0*diag(d)
        Backward.var.new[[t]][[i]] =  0*diag(d)
      }
    }
    
    
    ### update of sigma_beta
    
    sigma_beta_new = (sigma_beta_prior^(-2)+ V_beta_cumulative )^(-1/2)
    
    
    ### update of mean_beta
    
    mean_beta_new = as.numeric(sigma_beta_new^2* (sigma_beta_prior^(-2)*mean_beta_prior + M_beta_cumulative) )
    
    # mean_beta_new = 0
    
    ### update of Sigma_X, Mean_X
    
    

    for (i in 1:n){
      Mean_X_para =NULL
      Sigma_X_para = NULL
      for(t in 1:T){
        for (j in 1:n){
          if (j > i){
            M_i = Mean_X[[t]][i,]
            M_j = Mean_X[[t]][j,]
            V_i = Sigma_X[[t]][[i]]
            V_j = Sigma_X[[t]][[j]]
          } else if (j < i)
          {
            M_i = Mean_X[[t]][i,]
            M_j = Mean_X_new[[t]][j,]
            V_i = Sigma_X[[t]][[i]]
            V_j = Sigma_X_new[[t]][[j]]
          }
          
          if (j != i){
            V_X_cumulative[[t]][[i]] = V_X_cumulative[[t]][[i]] -2* Xi[[t]][i,j]*( M_j %*% t(M_j) + V_j)*alpha
            M_X_cumulative[[t]][i,] = M_X_cumulative[[t]][i,] +
              (Y[[t]][i,j]-0.5+2*Xi[[t]][i,j]*mean_beta_new) * M_j*alpha
          }
        }
      }
     
      
      Forward.var.new[[1]][[i]] = -rho^(-2)*trans_sd^(4)*(sigma_X1^(-2)*diag(d)+rho^(2)*trans_sd^(-2)*diag(d)+V_X_cumulative[[1]][[i]])
      
      Forward.Mean.new[[1]][i,] = -rho^(-1)*trans_sd^2*
        (M_X_cumulative[[1]][i,])
      
      
      
      Backward.var.new[[T]][[i]] = -rho^(-2)*trans_sd^(4)*(trans_sd^(-2)*diag(d)+V_X_cumulative[[T]][[i]])
      
      Backward.Mean.new[[T]][i,] = -rho^(-1)*trans_sd^2*
        (M_X_cumulative[[T]][i,])
      
      if (T>2){
        
        for (t2 in 2:(T-1)){
          
          
          
          
          Forward.var.new[[t2]][[i]] = -rho^(-2)*trans_sd^(4)*(ginv(Forward.var.new[[t2-1]][[i]])+(1+rho^2)*trans_sd^(-2)*diag(d)+V_X_cumulative[[t2]][[i]])
          
          Forward.Mean.new[[t2]][i,] = -rho^(-1)*trans_sd^2*
            (ginv(Forward.var.new[[t2-1]][[i]])%*%  Forward.Mean.new[[t2-1]][i,]  +M_X_cumulative[[t2]][i,])
          
          
          Backward.var.new[[T+1-t2]][[i]] = -rho^(-2)*trans_sd^(4)*(ginv(Backward.var.new[[T+2-t2]][[i]])+(1+rho^2)*trans_sd^(-2)*diag(d)+V_X_cumulative[[T+1-t2]][[i]])
          
          Backward.Mean.new[[T+1-t2]][i,] = -rho^(-1)*trans_sd^2*
            (ginv(Backward.var.new[[T+2-t2]][[i]])%*% Backward.Mean.new[[T+2-t2]][i,] +M_X_cumulative[[T+1-t2]][i,])
          
          
          
        }
      }

      
      for(t in 1:T){
        if(t == 1)
        { Sigma_X_new[[t]][[i]] = ginv(sigma_X1^(-2)*diag(d)+rho^2*trans_sd^(-2)*diag(d)+ginv(Backward.var.new[[t+1]][[i]])+V_X_cumulative[[t]][[i]])
        
        Sigma_X_new[[t]][[i]] = (Sigma_X_new[[t]][[i]] +t(Sigma_X_new[[t]][[i]]))/2
        
        Mean_X_new[[t]][i,] = Sigma_X_new[[t]][[i]]%*%
          (ginv(Backward.var.new[[t+1]][[i]])%*% Backward.Mean.new[[t+1]][i,]+M_X_cumulative[[t]][i,])
        }
        
        if(t==T)
        { Sigma_X_new[[t]][[i]] = ginv(trans_sd^(-2)*diag(d)+ginv(Forward.var.new[[t-1]][[i]])+V_X_cumulative[[t]][[i]])
        
        Sigma_X_new[[t]][[i]] = (Sigma_X_new[[t]][[i]] +t(Sigma_X_new[[t]][[i]]))/2
        
        Mean_X_new[[t]][i,] = Sigma_X_new[[t]][[i]]%*%
          (ginv(Forward.var.new[[t-1]][[i]])%*%  Forward.Mean.new[[t-1]][i,]+M_X_cumulative[[t]][i,])
        }
        
        if(1<t  && t<T)
        { Sigma_X_new[[t]][[i]] = ginv((1+rho^2)*trans_sd^(-2)*diag(d)+ginv(Forward.var.new[[t-1]][[i]])+ginv(Backward.var.new[[t+1]][[i]])+V_X_cumulative[[t]][[i]])
        
        Sigma_X_new[[t]][[i]] = (Sigma_X_new[[t]][[i]] +t(Sigma_X_new[[t]][[i]]))/2
        
        Mean_X_new[[t]][i,] = Sigma_X_new[[t]][[i]]%*%
          (ginv(Forward.var.new[[t-1]][[i]])%*%  Forward.Mean.new[[t-1]][i,]+
             ginv(Backward.var.new[[t+1]][[i]])%*% Backward.Mean.new[[t+1]][i,]+M_X_cumulative[[t]][i,])
        }
        
      }
      
      Cross_cov[[1]][[i]][1:d,1:d] = sigma_X1^(-2)*diag(d)+trans_sd^(-2)*diag(d)+V_X_cumulative[[1]][[i]]
      Cross_cov[[1]][[i]][1:d,(d+1):(2*d)] = -(trans_sd^(-2))*diag(d)
      Cross_cov[[1]][[i]][(d+1):(2*d),1:d] = -(trans_sd^(-2))*diag(d)
      Cross_cov[[1]][[i]][(d+1):(2*d),(d+1):(2*d)] = (trans_sd^(-2))*diag(d)+ginv(Backward.var.new[[3]][[i]])+(trans_sd^(-2))*diag(d)+V_X_cumulative[[2]][[i]]
      Cross_cov[[1]][[i]] = ginv(Cross_cov[[1]][[i]])
      
      for ( t in 2:(T-2)){
        Cross_cov[[t]][[i]][1:d,1:d] = (trans_sd^(-2))*diag(d)+ ginv(Forward.var.new[[t-1]][[i]]) +(trans_sd^(-2))*diag(d)+V_X_cumulative[[t]][[i]]
        Cross_cov[[t]][[i]][1:d,(d+1):(2*d)] = -(trans_sd^(-2))*diag(d)
        Cross_cov[[t]][[i]][(d+1):(2*d),1:d] = -(trans_sd^(-2))*diag(d)
        Cross_cov[[t]][[i]][(d+1):(2*d),(d+1):(2*d)] = (trans_sd^(-2))*diag(d)+ginv(Backward.var.new[[t+2]][[i]])+(trans_sd^(-2))*diag(d)+V_X_cumulative[[t+1]][[i]]
        Cross_cov[[t]][[i]] = ginv(Cross_cov[[t]][[i]])
      }
      
      Cross_cov[[T-1]][[i]][1:d,1:d] = (trans_sd^(-2))*diag(d)+ ginv(Forward.var.new[[T-2]][[i]]) +(trans_sd^(-2))*diag(d)+V_X_cumulative[[T-1]][[i]]
      Cross_cov[[T-1]][[i]][1:d,(d+1):(2*d)] = -(trans_sd^(-2))*diag(d)
      Cross_cov[[T-1]][[i]][(d+1):(2*d),1:d] = -(trans_sd^(-2))*diag(d)
      Cross_cov[[T-1]][[i]][(d+1):(2*d),(d+1):(2*d)] = (trans_sd^(-2))*diag(d)+V_X_cumulative[[T]][[i]]
      Cross_cov[[T-1]][[i]] = ginv(Cross_cov[[T-1]][[i]])
    }
    
    b=0
    b0 = 0
    for (i0 in 1:n){
      for (t in 1:(T-1)){
        
        EX = sum(Mean_X_new[[t+1]][i0,]^2)+sum(Mean_X_new[[t]][i0,]^2)-2*t(Mean_X_new[[t]][i0,]) %*% Mean_X_new[[t+1]][i0,]+
          sum(diag(x= as.matrix(Cross_cov[[t]][[i0]][1:d,1:d]+ Cross_cov[[t]][[i0]][(d+1):(2*d),(d+1):(2*d)] -
                                  2*Cross_cov[[t]][[i0]][1:d,(d+1):(2*d)])))
        
        # EX = sum(Mean_X_new[[t+1]][i,]^2)+sum(Mean_X_new[[t]][i,]^2)-2*t(Mean_X_new[[t]][i,]) %*% Mean_X_new[[t+1]][i,]+
        #   sum(diag(x=Sigma_X_new[[t]][[i]]))+  sum(diag(x = Sigma_X_new[[t+1]][[i]])) 
        
        b = b+ EX
      }
      b0 = b0 + sum(Mean_X_new[[1]][i0,]^2)+sum(diag(x=Sigma_X_new[[1]][[i0]]))
    }
    
    ded_free = n*(T-1)*d/2-1
    
    # if (abs(ded_free)<200){
    #   tau_inv_sq = besselK(x=sqrt(b), nu=ded_free+1)/sqrt(b)/besselK(x=sqrt(b), nu=ded_free)-2*ded_free/b 
    # }  else{
    #   tau_inv_sq =  exp(besselK.nuAsym(x=sqrt(b), nu=-ded_free-1,k.max=5,log=T)-besselK.nuAsym(x=sqrt(b), nu=-ded_free,k.max=5,log=T))/sqrt(b)-2*ded_free/b
    # }
    tau_inv_sq = (n*(T-1)*d+1)/(b+1)
    
    trans_sd = as.numeric(1/sqrt(tau_inv_sq))
    
    #  cat(trans_sd,'\n')
    
    sg_inve_sq = (n*d+1)/(b0+1)
    
    sigma_X1 = as.numeric(1/sqrt(sg_inve_sq))

    
    
    ### evaluation through roc 
    
    for (t in 1:T){
      auc_mean = rep((n-1)*n,0)
      auc_res = rep((n-1)*n,0)
      r=1
      for (i in 1:n){
        for (j in 1:n){
          if (j<i){
            auc_mean[r] = 1/(1+exp(-mean_beta_new-t( Mean_X_new[[t]][i,])%*% Mean_X_new[[t]][j,]))
            auc_res[r] = Y[[t]][i,j]
            r=r+1
          }
        }
      }
      roc_obj <- roc(auc_res, auc_mean,quiet=TRUE)
      
      train_auc[k] = train_auc[k]+auc(roc_obj)
    }
    
    train_auc[k] = train_auc[k]/T
    
    
    
 
    if(train_auc[k]-train_auc[k-1]< gap && k>min_iter+1){
      ind = 1
    }    else{
      cat(train_auc[k],'\n')

 
      
      Mean_X = Mean_X_new
      
      Sigma_X = Sigma_X_new
      
      mean_beta = as.numeric( mean_beta_new)
      
      sigma_beta = sigma_beta_new
      
      
      

      
      
      k = k+1
    }
    
  }
  
  
  return(list(norm=train_auc[2:(k-1)],trans_sd=trans_sd, Mean_X= Mean_X, Sigma_X= Sigma_X,  mean_beta  = mean_beta, sigma_beta =sigma_beta , iter = k-1 , AUC =train_auc[k-1]))
}