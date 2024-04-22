library(invgamma)


mcmc_DN_adaptive = function(Y,mean_beta_prior, d=2,sigma_beta_prior,  rho=1, gap =50,min_iter=1,alpha=0.95,max_iter=100){
  

  
  T = length(Y)
  
  n = length(Y[[1]][1,])
  
  
  #--------------------------------Model Initialization---------------------
  
  #### target parameters
  
  
  
  mean_beta = 0;    #mean of beta, 
  
  Mean_X = vector("list", T); #Mean of X, Mean_X[[t]][i,] 
  X_out = vector("list", T);
  
  for (t in 1:T){
    Mean_X[[t]] = matrix(rep(0,n*d),nrow = n)
    X_out[[t]] = matrix(rep(0,n*d),nrow = n)
    for (i in 1:n){
      Mean_X[[t]][i,] =  rmvnorm(n=1,mean=rep(0,d), sigma=diag(d)) # randomly initialization for mu_{it}
    }
  }

  
  #### tangent parameters
  
  Xi = vector("list", T)

  
  for (t in 1:T) { 
    Xi[[t]] = matrix(rep(1,n*n),nrow = n)
  }
  
  
  beta_out = 0
  #--------------------------------Algorithm---------------------
  K= 50000
  
  err =rep(0,K)

  train_auc =rep(0,K)
  
  train_auc[1]=-10
  
  k = 2
  
  ind = 0
  
  tau = rep(1,n)
  
  v0 = tau
  
  lambda = matrix(rep(0.5,n*(T-1)),ncol=T-1)
  
  v = lambda
  
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
          Xi[[t]][i,j] = rpg(num=1, h=1, z=t(M_i)%*%M_j +mean_beta)
          if (j!= i){
            V_beta_cumulative= V_beta_cumulative +  Xi[[t]][i,j]*alpha
            M_beta_cumulative = M_beta_cumulative + (Y[[t]][i,j]-0.5-Xi[[t]][i,j]* t(M_i) %*% M_j)*alpha
          }
        }
      }
    }
    
    mean_beta_new = 0;    
    
    sigma_beta_new = 0  
    
    Mean_X_new = vector("list", T)
    
    Sigma_X_new = replicate(n=T, expr=list())
    
    
    for (t in 1:T){
      Mean_X_new[[t]] = matrix(rep(0,n*d),nrow = n)
      Sigma_X_new[[t]] = vector("list", n)
      for (i in 1:n){
        Sigma_X_new[[t]][[i]] =  0*diag(d)         
      }
    }
    
    
    ### update of sigma_beta
    
    sigma_beta_new = (sigma_beta_prior^(-2)+ V_beta_cumulative )^(-1/2)
    
    
    ### update of mean_beta
    
    mean_beta_new = as.numeric(sigma_beta_new^2* (sigma_beta_prior^(-2)*mean_beta_prior + M_beta_cumulative) )
    
    
    mean_beta = rnorm(n=1, mean= mean_beta_new,sd = sigma_beta_new)
    
    
    ### update of Sigma_X, Mean_X
    
    for (t in 1:T){ tryCatch({
      for (i in 1:n){
        for (j in 1:n){
          if (j > i){
            M_i = Mean_X[[t]][i,]
            M_j = Mean_X[[t]][j,]
          } else if (j < i)
          {
            M_i = Mean_X[[t]][i,]
            M_j = Mean_X[[t]][j,]
          }
          
          if (j!= i){
            V_X_cumulative[[t]][[i]] = V_X_cumulative[[t]][[i]] + Xi[[t]][i,j]*( M_j %*% t(M_j) )*alpha
            M_X_cumulative[[t]][i,] = M_X_cumulative[[t]][i,] +
              (Y[[t]][i,j]-0.5-Xi[[t]][i,j]*mean_beta) * M_j*alpha
          }
        }
        
        if(t == 1)
        { Sigma_X_new[[t]][[i]] = ginv(sigma_X1^(-2)*diag(d)+(lambda[i,t]*tau[i])^(-1)*diag(d)+V_X_cumulative[[t]][[i]])
        
        Sigma_X_new[[t]][[i]] = (Sigma_X_new[[t]][[i]] +t(Sigma_X_new[[t]][[i]]))/2
        
        Mean_X_new[[t]][i,] = Sigma_X_new[[t]][[i]]%*%
          ((lambda[i,t]*tau[i])^(-1)*Mean_X[[t+1]][i,]+M_X_cumulative[[t]][i,])
        }
        
        
        if(1<t  && t<T)
        { Sigma_X_new[[t]][[i]] = ginv((lambda[i,t-1]*tau[i])^(-1)*diag(d)+(lambda[i,t]*tau[i])^(-1)*diag(d)+V_X_cumulative[[t]][[i]])
        
        Sigma_X_new[[t]][[i]] = (Sigma_X_new[[t]][[i]] +t(Sigma_X_new[[t]][[i]]))/2
        
        Mean_X_new[[t]][i,] = Sigma_X_new[[t]][[i]]%*%
          (rho*(lambda[i,t-1]*tau[i])^(-1)*Mean_X[[t-1]][i,]+(lambda[i,t]*tau[i])^(-1)*Mean_X[[t+1]][i,]+M_X_cumulative[[t]][i,])
        }
        
        if(t==T)
        { Sigma_X_new[[t]][[i]] = ginv((lambda[i,t-1]*tau[i])^(-1)*diag(d)+V_X_cumulative[[t]][[i]])
        
        Sigma_X_new[[t]][[i]] = (Sigma_X_new[[t]][[i]] +t(Sigma_X_new[[t]][[i]]))/2
        
        Mean_X_new[[t]][i,] = Sigma_X_new[[t]][[i]]%*%
          ((lambda[i,t-1]*tau[i])^(-1)*diag(d)*Mean_X[[t-1]][i,]+M_X_cumulative[[t]][i,])
        }
        
        Mean_X[[t]][i,] =rmvnorm(n=1, mean = Mean_X_new[[t]][i,], sigma = Sigma_X_new[[t]][[i]])
      }
      


    },error=function(e){})
    }
    
    for (t in 1:(T-1)){
      for (i in 1:n){
      v[i,t] = rinvgamma(n=1,shape=1,rate = 1+1/lambda[i,t])
      }
    }
    EX=rep(0,t-1)
    b0 = 0 
    for (i in 1:n){
    b = 0
    EX=rep(0,t-1)
    for (t in 1:(T-1)){
      EX[t] = sum(Mean_X[[t+1]][i,]^2)+sum(Mean_X[[t]][i,]^2)-2*t(Mean_X[[t]][i,]) %*% Mean_X[[t+1]][i,]
      lambda[i,t] =rinvgamma(n=1, shape=(d+1)/2,  rate=1/v[i,t] + EX[t]/(2*tau[i]))
      b = b+ EX[t]/lambda[i,t]
    }
    b0 = b0 + sum(Mean_X[[1]][i,]^2)
    
    
    
    v0 = rinvgamma(n=1,shape=1,rate = 1+1/tau)
    
    tau[i]  =rinvgamma(n=1, shape=d*(T-1)/2+1/2,  rate=1/v0+b/2) 
    }
    sg_inve_sq = rinvgamma(n=1,shape= (n*d+1)/2, rate=(b0+1)/2)
    
    sigma_X1 = as.numeric(sqrt(sg_inve_sq))

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
      roc_obj <- roc(response =auc_res, predictor =auc_mean,quiet=TRUE)
      
      train_auc[k] = train_auc[k]+auc(roc_obj)
    }
    train_auc[k] = train_auc[k]/T
    
    if(  k>max_iter){
      #(train_auc[k]-train_auc[k-1]< 0.01)
      ind = 1

    }    else{
    
      
      cat(train_auc[k],'\n')
      
      cat(k-1,'\n')

      
      if (k > (max_iter/2)){
        for(t in 1:T){
          X_out[[t]] =  X_out[[t]]+Mean_X[[t]]/(max_iter/2)
        }
        beta_out = beta_out + mean_beta/(max_iter/2)
      }
    

      
      
      k = k+1
    
    }
    
  
  }
  
  return(list(norm=train_auc[2:(k-1)], Mean_X= X_out, mean_beta  = beta_out,  iter = k-1 , AUC =train_auc[k-1]))
}
