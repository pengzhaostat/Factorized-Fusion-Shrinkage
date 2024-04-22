#----------------------------------------------------------------------------
# Gaussian matrix factorization via SMF variational inference
# Input Variables:

# Y: Observed values, list of length T, Y[[t]]: n*p matrix values at time point t
# mean_beta_prior,sigma_beta_prior: mean and sd for prior of beta
# gap_per_iter: gap between errors for convergence
# max_iter: maximal cycles in computation
# alpha: fractional power of the likelihood, fixed to be 0.95 in this paper
# global_prior: types of global_prior, Cauthy, fixed or Gamma
# type: types of observations, full, off-diagonal

#Output variables:
# err: dynamic of the training RMSE
# Mean_X: variational mean of each subject of X, list of length T, Mean_X[[t]]: n*d matrix at time point t
# Sigma_X: variational covariance matrix of each subject of X
# Mean_Z: variational mean of each subject of Z
# Sigma_Z: variational covariance matrix of each subject of Z
# mean_beta: variational mean of intercept beta
# sigma_beta: variational sd of intercept beta 
# iter: k-1
# tau_1: means of global_prior of X
# tau_2: means of global_prior of Z
# lambda_X: means of local prior of X
# lambda_Z: means of local prior of Z


MP_Gaussain_weighted_adaptive = function(Y,tau=0.01, sigma=0.5 , gap =1e-6, max_iter=20, X_init = NULL, 
                                alpha=0.95, d=2,mean_beta_prior=0, sigma_beta_prior=10, type = 'full',
                                global_prior=global_prior,gap_per_iter=1e-3){
  
  tau_1 = 1/tau^2
  
  tau_2 = 1/tau^2
  
  T = length(Y)
  
  n = length(Y[[1]][,1])
  
  p = length(Y[[1]][1,])
  
  #--------------------------------Model Initialization---------------------
  mean_beta = 0;    #mean of beta, 
  
  sigma_beta = 1;  #sd of beta, 
  
  #### target parameters
  Mean_X = vector("list", n); #Mean of X, Mean_X[[t]][i,] 
  
  Mean_Z = vector("list", p);
  
  Sigma_X = replicate(n=n, expr=list()); # covariance matrix of X, Sigma_x[[t]][[i]] is the covariance matrix Sigma_{it}
  
  for (i in 1:n){
    Mean_X[[i]] = matrix(rnorm(d*T),nrow = T)
    Sigma_X[[i]] = vector("list", T)
    for (t in 1:T){
      Sigma_X[[i]][[t]] =  0.1*diag(d)          # initialization for Sigma_{it}
    }
  }
  
  Sigma_Z = replicate(n=p, expr=list())
  
  
  for (i in 1:p){
    Mean_Z[[i]] = matrix(rep(0.1,d*T),nrow = T)
    Sigma_Z[[i]] = vector("list", T)
    for (t in 1:T){
      Sigma_Z[[i]][[t]] =  0.1*diag(d)          # initialization for Sigma_{it}
    }
  }
  
  Mean_X_new = vector("list", n)
  
  Sigma_X_new = replicate(n=n, expr=list())
  
  Mean_Z_new = vector("list", p)
  
  Sigma_Z_new = replicate(n=p, expr=list())
  
  lambda_X = matrix(rep(1,n*(T-1)),nrow=n)

  v_X = matrix(rep(1,n*(T-1)),nrow=n)
  
  lambda_Z = matrix(rep(1,p*(T-1)),nrow=p)
  
  v_Z = matrix(rep(1,p*(T-1)),nrow=p)
  
  EX = rep(0,T-1)
  
  Mean_X_new = Mean_X
  for (i in 1:n){
    Sigma_X_new[[i]] = vector("list", T)
    for (t in 1:T){
      Sigma_X_new[[i]][[t]] =  diag(d)
    }
  }
  
  Mean_Z_new = Mean_Z
  for (i in 1:p){
    Sigma_Z_new[[i]] = vector("list", T)
    for (t in 1:T){
      Sigma_Z_new[[i]][[t]] =  diag(d)
    }
  }
  
  sigma_0_X = rep(0.5,n)
  
  sigma_0_Z =rep(0.5,p)
  
  #--------------------------------Algorithm---------------------
  K= 1000
  
  err =rep(0,K)
  
  err[1] =1e10
  
  norm_stop = rep(0,K)
  
  k = 2
  
  ind = 0
  
  err_node =1:n
  
  


  if (!is.null(X_init)){
    Mean_X = X_init
  }
  
  while(k<K && ind ==0){
    
    V_beta_cumulative = 0
     
    M_beta_cumulative = 0
    
    V_X_cumulative =  replicate(n=n, expr=list())
    
    M_X_cumulative = vector("list", n)
    
    
    for (i in 1:n){
      M_X_cumulative[[i]] = matrix(rep(0,T*d),nrow = T)
      V_X_cumulative[[i]] = vector("list", T)
      for (t in 1:T){
        V_X_cumulative[[i]][[t]] =  0*diag(d)        
      }
    }
    
    for (t in 1:T){
      for( i in 1:n){
        for (j in 1:p){
          M_i = Mean_X_new[[i]][t,]
          M_j = Mean_Z_new[[j]][t,]
          V_i = Sigma_X_new[[i]][[t]]
          V_j = Sigma_Z_new[[j]][[t]]
            M_beta_cumulative = M_beta_cumulative + (Y[[t]][i,j]- t(M_i) %*% M_j)*alpha
        }
      }
    }
    ### update of sigma_beta
    
  #  sigma_beta_new = (sigma_beta_prior^(-2)+ T*n*p*sigma^(-2) )^(-1/2)
    
    sigma_beta_new =1
    
    ### update of mean_beta
    
  #  mean_beta_new = as.numeric(sigma_beta_new^2* (sigma_beta_prior^(-2)*mean_beta_prior + sigma^(-2)*M_beta_cumulative) )
    
    mean_beta_new = 0 
    

    ### update of Sigma_X, Mean_X
    
    if (type == 'full'){
    
    for (i in 1:n){
      for (t in 1:T){
        for (j in 1:p){
            M_j = Mean_Z_new[[j]][t,]
            V_j = Sigma_Z_new[[j]][[t]]
            V_X_cumulative[[i]][[t]] = V_X_cumulative[[i]][[t]] + sigma^(-2)*(M_j %*% t(M_j) + V_j) *alpha
            M_X_cumulative[[i]][t,] = M_X_cumulative[[i]][t,] +   sigma^(-2)*(Y[[t]][i,j]-mean_beta_new) * M_j *alpha
        }
        M_X_cumulative[[i]][t,] = solve(V_X_cumulative[[i]][[t]]) %*% M_X_cumulative[[i]][t,]
        }
        
       MF_gibbs =  MP_gibbs_mult_Sigma (Y = M_X_cumulative[[i]], init_sigma_X1=sigma_0_X[i], V_X_cumulative[[i]],  tau =tau_1, max_iter = 100, 
                                         gap = gap_per_iter,global_prior=global_prior)

       Mean_X_new[[i]] = MF_gibbs$Mean
       Sigma_X_new[[i]] = MF_gibbs$Sigma
       lambda_X[i,] = MF_gibbs$lambda
       sigma_0_X[i] = MF_gibbs$sigma_X1
    }
    } else if (type == 'upper'){
      
      for (i in 1:n){
        for (t in 1:T){
          for (j in 1:p){
            if (j>=i){
            M_j = Mean_Z_new[[j]][t,]
            V_j = Sigma_Z_new[[j]][[t]]
            V_X_cumulative[[i]][[t]] = V_X_cumulative[[i]][[t]] + sigma^(-2)*(M_j %*% t(M_j) + V_j) *alpha
            M_X_cumulative[[i]][t,] = M_X_cumulative[[i]][t,] +   sigma^(-2)*(Y[[t]][i,j]-mean_beta_new) * M_j *alpha
            }
          }
          M_X_cumulative[[i]][t,] = solve(V_X_cumulative[[i]][[t]]) %*% M_X_cumulative[[i]][t,]
        }
        
        MF_gibbs =  MP_gibbs_mult_Sigma (Y = M_X_cumulative[[i]], init_sigma_X1=sigma_0_X[i], V_X_cumulative[[i]],  tau =tau_1, max_iter = 100, 
                                          gap = gap_per_iter,global_prior=global_prior)
        
        Mean_X_new[[i]] = MF_gibbs$Mean
        Sigma_X_new[[i]] = MF_gibbs$Sigma
        lambda_X[i,] = MF_gibbs$lambda
        sigma_0_X[i] = MF_gibbs$sigma_X1
      }
    }
    
    
    
    V_Z_cumulative =  replicate(n=p, expr=list())
    
    M_Z_cumulative = vector("list", p)
    
    
    for (i in 1:p){
      M_Z_cumulative[[i]] = matrix(rep(0,T*d),nrow = T)
      V_Z_cumulative[[i]] = vector("list", T)
      for (t in 1:T){
        V_Z_cumulative[[i]][[t]] =  0*diag(d)        
      }
    }
    
    if (type == 'full'){
 
    for (j in 1:p){
      for (t in 1:T){
        for (i in 1:n){
          M_i = Mean_X_new[[i]][t,]
          V_i = Sigma_X_new[[i]][[t]]
          V_Z_cumulative[[j]][[t]] = V_Z_cumulative[[j]][[t]] + sigma^(-2)*(M_i %*% t(M_i) + V_i) *alpha
          M_Z_cumulative[[j]][t,] = M_Z_cumulative[[j]][t,] +   sigma^(-2)*(Y[[t]][i,j]-mean_beta_new) * M_i *alpha
        }
        M_Z_cumulative[[j]][t,] = solve(V_Z_cumulative[[j]][[t]]) %*% M_Z_cumulative[[j]][t,]
      }
      
      MF_gibbs =  MP_gibbs_mult_Sigma (Y = M_Z_cumulative[[j]], init_sigma_X1=sigma_0_Z[j], V_Z_cumulative[[j]],  tau =tau_2, max_iter = 100, 
                                        gap = gap_per_iter,global_prior=global_prior)
      
      Mean_Z_new[[j]] = MF_gibbs$Mean
      Sigma_Z_new[[j]] = MF_gibbs$Sigma
      lambda_Z[j,] = MF_gibbs$lambda
      sigma_0_Z[j] = MF_gibbs$sigma_X1
    }
    } else if (type =='upper'){
      for (j in 1:p){
        for (t in 1:T){
          for (i in 1:n){
            if (j >=i){
            M_i = Mean_X_new[[i]][t,]
            V_i = Sigma_X_new[[i]][[t]]
            V_Z_cumulative[[j]][[t]] = V_Z_cumulative[[j]][[t]] + sigma^(-2)*(M_i %*% t(M_i) + V_i) *alpha
            M_Z_cumulative[[j]][t,] = M_Z_cumulative[[j]][t,] +   sigma^(-2)*(Y[[t]][i,j]-mean_beta_new) * M_i *alpha
            }
          }
          M_Z_cumulative[[j]][t,] = solve(V_Z_cumulative[[j]][[t]]) %*% M_Z_cumulative[[j]][t,]
        }
        
        MF_gibbs =  MP_gibbs_mult_Sigma (Y = M_Z_cumulative[[j]], init_sigma_X1=sigma_0_Z[j], V_Z_cumulative[[j]],  tau =tau_2, max_iter = 100, 
                                          gap = gap_per_iter,global_prior=global_prior)
        
        Mean_Z_new[[j]] = MF_gibbs$Mean
        Sigma_Z_new[[j]] = MF_gibbs$Sigma
        lambda_Z[j,] = MF_gibbs$lambda
        sigma_0_Z[j] = MF_gibbs$sigma_X1
      }
    }
    
    
    distance_squared_inner_prod = function(mu1,mu2,Sigma1,Sigma2){
      value = t(mu1)%*% mu2 %*% t(mu2) %*%mu1 +sum(diag(Sigma1 %*% Sigma2)) +  t(mu1)%*% Sigma2 %*%mu1 + t(mu2)%*% Sigma1 %*%mu2
      
      return(value)
    }
    
    if (type == 'full'){
    
    pred_mean = rep(T*p*n,0)
    res = rep(T*p*n,0)
    sigma_sum = 0
    r=1
 
      for (i in 1:n){
        for (j in 1:p){
            for (t in 1:T){
            pred_mean[r] = mean_beta_new+t(Mean_X_new[[i]][t,])%*% Mean_Z_new[[j]][t,]
            res[r] = Y[[t]][i,j]
            r=r+1
            
            sigma_sum =  sigma_sum + Y[[t]][i,j]^2-2*Y[[t]][i,j]*t(Mean_X_new[[i]][t,])%*% Mean_Z_new[[j]][t,] +
              distance_squared_inner_prod(Mean_X_new[[i]][t,],Mean_Z_new[[j]][t,],Sigma_X_new[[i]][[t]],Sigma_Z_new[[j]][[t]])
        }
      }
      }
    sg_inve_sq = (n*p*T+1)/(sigma_sum+1)
    
    sigma_new = as.numeric(1/sqrt(sg_inve_sq))
    } else if (type == 'upper'){
      
      pred_mean = rep(T*p*n,0)
      res = rep(T*p*n,0)
      r=1
      
      sigma_sum = 0
      for (i in 1:n){
        for (j in 1:p){
          for (t in 1:T){
            if (p>=n){
            pred_mean[r] = mean_beta_new+t(Mean_X_new[[i]][t,])%*% Mean_Z_new[[j]][t,]
            res[r] = Y[[t]][i,j]
            r=r+1
            
            sigma_sum =  sigma_sum + Y[[t]][i,j]^2-2*Y[[t]][i,j]*t(Mean_X_new[[i]][t,])%*% Mean_Z_new[[j]][t,] +
              distance_squared_inner_prod(Mean_X_new[[i]][t,],Mean_Z_new[[j]][t,],Sigma_X_new[[i]][[t]],Sigma_Z_new[[j]][[t]])
            }
          }
        }
      }
      
      sg_inve_sq = (n*p*T+1)/(sigma_sum+1)
      
      sigma_new = as.numeric(1/sqrt(sg_inve_sq))
    }
    
    
    

    
    
    err[k]=sqrt(sum((res-pred_mean)^2/T/n/p))
    
    
   
    cat(err[k],'\n')
    
    cat('iteration:',k-1,'\n')
  
    
    if(abs(err[k]-err[k-1])< gap ||  k>max_iter+1){
      ind = 1
    }    else{
      

      
      Mean_X = Mean_X_new
      
      Sigma_X = Sigma_X_new
      
      Mean_Z = Mean_Z_new
      
      Sigma_Z = Sigma_Z_new
      
      sigma_beta = sigma_beta_new
      
      mean_beta = mean_beta_new
      
      sigma =sigma_new
      
      
      k = k+1
    }
    
  }
  
  return(list(err=err[2:(k-1)], Mean_X= Mean_X, Sigma_X= Sigma_X, Mean_Z= Mean_Z, Sigma_Z= Sigma_Z,  iter = k-1, sigma= sigma,
              tau_1=tau_1, tau_2 =tau_2,lambda_X =lambda_X,lambda_Z = lambda_Z, mean_beta=mean_beta,sigma_beta=sigma_beta, sigma_0_X = sigma_0_X, sigma_0_Z=sigma_0_Z))
}
