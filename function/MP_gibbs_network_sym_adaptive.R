distance_squared_inner_prod = function(mu1,mu2,Sigma1,Sigma2){
  value = t(mu1)%*% mu2 %*% t(mu2) %*%mu1 +sum(diag(Sigma1 %*% Sigma2)) +  t(mu1)%*% Sigma2 %*%mu1 + t(mu2)%*% Sigma1 %*%mu2
  
  return(value)
}

# Core function to perform group_wise fused shirnkage for a single subject, all other functions, like Gaussian matrix factorization, binary network and tensor models will use this function iteratively

MP_binary_weighted_adaptive = function(Y,tau=0.01, gap =0.01, max_iter=2000, X_init = NULL,alpha = 0.95, d=2,mean_beta_prior=0,
                                       sigma_beta_prior=sqrt(10),global_prior='Cauthy'){
  
  tau = 1/tau^2
  
  T = length(Y)
  
  n = length(Y[[1]][,1])
  
  d = 2
  #--------------------------------Model Initialization---------------------
  mean_beta = 0;    #mean of beta, 
  
  sigma_beta = 1;  #sd of beta, 
  
  #### target parameters
  Mean_X = vector("list", n); #Mean of X, Mean_X[[t]][i,] 
  
  Sigma_X = replicate(n=n, expr=list()); # covariance matrix of X, Sigma_x[[t]][[i]] is the covariance matrix Sigma_{it}
  
  for (i in 1:n){
    Mean_X[[i]] = 0.1*matrix(rnorm(T*d),nrow = T)
    Sigma_X[[i]] = vector("list", T)
    for (t in 1:T){
      Sigma_X[[i]][[t]] =  0*diag(d)          # initialization for Sigma_{it}
    }
  }
  
 
  Mean_X_new = vector("list", n)
  
  Sigma_X_new = replicate(n=n, expr=list())
  

  lambda_X = matrix(rep(1,n*(T-1)),nrow=n)
  
  v_X = matrix(rep(1,n*(T-1)),nrow=n)
  

  
  EX = rep(0,T-1)
  
  for (i in 1:n){
    Mean_X_new[[i]] = matrix(rep(0,T*d),nrow = T)
    Sigma_X_new[[i]] = vector("list", T)
    for (t in 1:T){
      Sigma_X_new[[i]][[t]] =  0*diag(d)
    }
  }

  
  sigma_0_X = rep(0.5,n)
  

  
  Xi = vector("list", T)

  
  for (t in 1:T) { 
    Xi[[t]] = matrix(rep(1,n*n),nrow = n)
  }
  #--------------------------------Algorithm---------------------
  K= 1000
  
  err =rep(0,K)
  
  err[1] =0
  
  norm_stop = rep(0,K)
  
  k = 2
  
  ind = 0
  
  
  EX = rep(0,T-1)
  
  if (!is.null(X_init)){
    Mean_X = X_init
  }
  
  
  while(k<K && ind ==0){
    
     V_beta_cumulative = 0
    # 
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
    
    # 
    for (t in 1:T){
      for( i in 1:n){
        for (j in 1:n){
          M_i = Mean_X[[i]][t,]
          M_j = Mean_X[[j]][t,]
          V_i = Sigma_X[[i]][[t]]
          V_j = Sigma_X[[j]][[t]]
          c = t(M_i) %*% M_j+ mean_beta
          Xi[[t]][i,j] = 1/(2*c)*(exp(c)-1)/(exp(c)+1)/(-2)
          if (j!= i){
            V_beta_cumulative= V_beta_cumulative -  2*Xi[[t]][i,j]*alpha
            M_beta_cumulative = M_beta_cumulative + (Y[[t]][i,j]-0.5+2*Xi[[t]][i,j]* t(M_i) %*% M_j)*alpha
          }
        }
      }
    }
    ### update of sigma_beta
    
 #   sigma_beta_new = (sigma_beta_prior^(-2)+ V_beta_cumulative )^(-1/2)
    
    sigma_beta_new = 0
    ### update of mean_beta
    
  #   mean_beta_new = as.numeric(sigma_beta_new^2* (sigma_beta_prior^(-2)*mean_beta_prior + M_beta_cumulative) )

    mean_beta_new = 0
    
    ### update of Sigma_X, Mean_X
    
    for (i in 1:n){
      for (t in 1:T){
        for (j in 1:n){
          if (j > i){
            M_i = Mean_X[[i]][t,]
            M_j = Mean_X[[j]][t,]
            V_i = Sigma_X[[i]][[t]]
            V_j = Sigma_X[[j]][[t]]
            V_X_cumulative[[i]][[t]] = V_X_cumulative[[i]][[t]] -2* Xi[[t]][i,j]*( M_j %*% t(M_j) + V_j)*alpha
            M_X_cumulative[[i]][t,] = M_X_cumulative[[i]][t,] +  (Y[[t]][i,j]-0.5+2*Xi[[t]][i,j]*mean_beta_new) * M_j*alpha
          } else if (j < i)
          {
            M_i = Mean_X[[i]][t,]
            M_j = Mean_X_new[[j]][t,]
            V_i = Sigma_X[[i]][[t]]
            V_j = Sigma_X_new[[j]][[t]] 
            V_X_cumulative[[i]][[t]] = V_X_cumulative[[i]][[t]] -2* Xi[[t]][i,j]*( M_j %*% t(M_j) + V_j)*alpha
            M_X_cumulative[[i]][t,] = M_X_cumulative[[i]][t,] + (Y[[t]][i,j]-0.5+2*Xi[[t]][i,j]*mean_beta_new)  * M_j*alpha
          }
        }
        M_X_cumulative[[i]][t,] = ginv(V_X_cumulative[[i]][[t]]) %*% M_X_cumulative[[i]][t,]
        }
        
      # MF_gibbs =  MP_gibbs_mult_Sigma (Y = M_X_cumulative[[i]], init_sigma_X1=sigma_0_X[i], V_X_cumulative[[i]],  tau =tau, max_iter = 100, 
      #                                  init_Mean=Mean_X[[i]], init_Sigma = Sigma_X[[i]],  gap = 1e-6)
      
      MF_gibbs =  MP_gibbs_mult_Sigma (Y = M_X_cumulative[[i]], Sigma_list = V_X_cumulative[[i]], tau=tau,  max_iter = 100, 
                                         gap = 1e-4, global_prior=global_prior)
      
#      MF_m = MF_gibbs$Mean
      
#      MF_m [MF_m>10] =10
#      MF_m [MF_m<-10] = -10
      
      Mean_X_new[[i]] = MF_gibbs$Mean
      
      Sigma_X_new[[i]] = MF_gibbs$Sigma
      lambda_X[i,] = MF_gibbs$lambda
      sigma_0_X[i] = MF_gibbs$sigma_X1
      
    }
 

    
    pred_mean = rep(T*n*n,0)
    res = rep(T*n*n,0)
    r=1
 
      for (i in 1:n){
        for (j in 1:n){
            for (t in 1:T){
            pred_mean[r] = 1/(1+exp(-mean_beta_new-t(-Mean_X_new[[i]][t,])%*% Mean_X_new[[j]][t,]))
            res[r] = Y[[t]][i,j]
            r=r+1
        }
      }
      }
    
    roc_obj <- roc(res, pred_mean,quiet=TRUE)
    err[k]= auc(roc_obj)
     

    
    if( abs(err[k]-err[k-1] < gap)  || k>max_iter+1){
      ind = 1
    }    else{
      
      cat(err[k],'\n')
      cat('iteration:',k-1,'\n')
      
      Mean_X = Mean_X_new
      
      Sigma_X = Sigma_X_new
      
      sigma_beta = sigma_beta_new
      
      mean_beta = mean_beta_new
      
      
      k = k+1
    }
    
  }
  tau = MF_gibbs$tau
  return(list(norm = norm_stop[2:(k-1)], Mean_X= Mean_X, Sigma_X= Sigma_X, iter = k-1, 
              tau=tau,lambda_X =lambda_X, mean_beta=mean_beta,sigma_beta=sigma_beta, sigma_0_X = sigma_0_X,Xi=Xi))
}
