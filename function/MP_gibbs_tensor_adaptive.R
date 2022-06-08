distance_squared_inner_prod = function(mu1,mu2,Sigma1,Sigma2){
  value = t(mu1)%*% mu2 %*% t(mu2) %*%mu1 +sum(diag(Sigma1 %*% Sigma2)) +  t(mu1)%*% Sigma2 %*%mu1 + t(mu2)%*% Sigma1 %*%mu2
  
  return(value)
}


distance_squared_element_prod = function(X_Sigma1,mu2,Sigma2){
  value = diag(mu2) %*% X_Sigma1 %*% diag(mu2) + X_Sigma1*Sigma2
  return(value)
}

MP_Gaussain_tensor_adaptive = function(Y,tau=0.01, sigma=0.5 , gap =1e-6, max_iter=20, X_init = NULL, 
                                         alpha=0.95, d=2,
                                         global_prior='Cauthy',gap_per_iter=1e-3){
  
  
  T = dim(Y)[length(dim(Y))]
  
  n_m = dim(Y)[-length(dim(Y))]
  
  n = max(n_m)
  
  m = length(n_m)
  
  Y_tensor = as.tensor(Y)
  

  #--------------------------------Model Initialization---------------------

  #### target parameters
  Mean_X = 0.1*rand_tensor(c(n,d,m,T))@data #Mean of X
  
  Sigma_X = array(rep(c(1,0,0,1),),dim=c(d,d,n,m,T)); # covariance matrix of X
  
  
  Mean_X_new = Mean_X
  
  Sigma_X_new =  Sigma_X
  
  lambda_X = array(rep(1,),dim=c(n,m,T-1))
  
  v_X = array(rep(1,),dim=c(n,m,T-1))
  
  tau = array(rep(1,),dim=c(n,m))
  
  sigma_0_X = array(rep(0.5,),dim=c(n,m))
  
  #--------------------------------Algorithm---------------------
  K= 1000
  
  err =rep(0,K)
  
  err[1] =1e10
  
  norm_stop = rep(0,K)
  
  k = 2
  
  ind = 0
  
  tau = 1/tau^2

  
  while(k<K && ind ==0){
    
    cat(k,'\n')
    
    for (m0 in 1:m){
      
      if(m0 == 1){m1 = 2; m2 = 3}
      if(m0 == 2){m1 = 1; m2 = 3}
      if(m0 == 3){m1 = 1; m2 = 2}
      
    
    V_X_cumulative =  replicate(n=n_m[m0], expr=list())
    
    M_X_cumulative = vector("list", n_m[m0])
    
    
    for (i in 1:n_m[m0]){
      M_X_cumulative[[i]] = matrix(rep(0,T*d),nrow = T)
      V_X_cumulative[[i]] = vector("list", T)
      for (t in 1:T){
        V_X_cumulative[[i]][[t]] =  0*diag(d)        
      }
    }
    
    ### update of Sigma_X, Mean_X
    
      
      for (i in 1:n_m[m0]){
        for (t in 1:T){
          for (j1 in 1:n_m[m1]){
            M_j1 = Mean_X_new[j1,,m1,t]
            V_j1 = Sigma_X_new[,,j1,m1,t]
            for (j2 in 1:n_m[m2]){
            M_j2 = Mean_X_new[j2,,m2,t]
            V_j2 = Sigma_X_new[,,j2,m2,t]
            V_X_cumulative[[i]][[t]] = V_X_cumulative[[i]][[t]] + sigma^(-2)*
              (distance_squared_element_prod(M_j1%*%t(M_j1)+V_j1,M_j2,V_j2)) *alpha
            M_X_cumulative[[i]][t,] = M_X_cumulative[[i]][t,] +   sigma^(-2)*
               (cs_unfold(Y_tensor[,,,t],m0)@data[(j2-1)*n_m[m1]+j1,i])* (M_j1*M_j2) *alpha
            }
          }
          M_X_cumulative[[i]][t,] = solve(V_X_cumulative[[i]][[t]]) %*% M_X_cumulative[[i]][t,]
          }
        
        cat(i,' ')
        
        MF_gibbs =  MP_gibbs_mult_Sigma (Y = M_X_cumulative[[i]],Sigma_list= V_X_cumulative[[i]], max_iter = 100, 
                                         gap = 1e-3)
  
     
        Mean_X_new[i,,m0,] = t(MF_gibbs$Mean)
        Sigma_X_new[,,i,m0,] = array(unlist(MF_gibbs$Sigma),dim=c(d,d,T))
        lambda_X[i,m0,] = MF_gibbs$lambda
        sigma_0_X[i,m0] = MF_gibbs$sigma_X1
      }
    }
      
    pred_mean = array(rep(0,),dim=c(n_m,T))
    for ( t in 1:T){
      for (d0 in 1:d){
        a = Mean_X_new[1:n_m[1],d0,1,t]
        for (i in 2:m){
          a = outer(a,Mean_X_new[1:n_m[i],d0,i,t])
        }
        pred_mean[,,,t] =  pred_mean[,,,t] +a
      }
    }
    
    
 
  
    
    err[k]=sqrt(sum((as.vector(pred_mean-Y))^2)/T/(n_m[1]*n_m[2]*n_m[3]))

    cat(err[k],'\n')
    
    cat('iteration:',k-1,'\n')
    
    
    if(abs(err[k]-err[k-1])< gap ||  k>max_iter+1){
      ind = 1
    }    else{
      
      Mean_X = Mean_X_new
      
      Sigma_X = Sigma_X_new
    
      
      
      k = k+1
    }
    
  }
  
  return(list(err=err[2:(k-1)], Mean_X= Mean_X, Sigma_X= Sigma_X,  iter = k-1, sigma= sigma,
              tau =tau,lambda_X =lambda_X, sigma_0_X = sigma_0_X,pred_mean=pred_mean))
}
