library(Bessel)

MP_gibbs_mult_Sigma =  function(Y, init_sigma_X1=NULL, Sigma_list, tau =1000, max_iter, init_Mean = NULL, init_Sigma = NULL,
                                init_lambda = NULL,gap =1e-6, global_prior='Cauthy'){
  
  T = length(Y[,1])
  d = length(Y[1,])
  K=2000
  X_Mean = matrix(rnorm(T*d),nrow = T)
  Forward.Mean.new = matrix(rep(0,T*d),nrow = T)
  Backward.Mean.new = matrix(rep(0,T*d),nrow = T)
  X_sigma = vector("list", T)
  Forward.var.new = vector("list", T)
  Backward.var.new = vector("list", T)
  Cross_cov = vector("list", T)
  mean_matrix = Y
  for (t in 1:T){
    X_sigma[[t]] = diag(d)
    Forward.var.new[[t]] = diag(d)
    Backward.var.new[[t]] = diag(d)
    Cross_cov [[t]] = diag(2*d)
    mean_matrix[t,] = Sigma_list[[t]] %*% Y[t,]
  }
  lambda = rep(1,T-1)
  lambda_new = rep(1,T-1)
  v = rep(1,T-1)
  
  tau_0 = 10000
  
  EX = rep(0,T-1)
  err =rep(0,K)
  
  err[1] =1e10
  v0 = 1
  k=2
  ind = 0
  
  sigma_X1 = 0.5
  
  X_Mean_old = X_Mean
  
  X_Sigma_old = X_sigma
  
  if (!is.null(init_Mean)){
    X_Mean = init_Mean
  }
  
  if (!is.null(init_Sigma)){
    X_sigma = init_Sigma
  }
  if (!is.null(init_lambda)){
    lambda = init_lambda
  }
  
  if (!is.null(init_sigma_X1)){
    sigma_X1 = init_sigma_X1
  }
  
  while(k<K && ind ==0){ 
    res <- try({
    
    Forward.var.new[[1]] = -(lambda[1]*tau)^(-2)*(sigma_X1^(-2)*diag(d)+(lambda[1]*tau)*diag(d)+Sigma_list[[1]])
    
    Forward.Mean.new[1,] = -(lambda[1]*tau)^(-1)*mean_matrix[1,]
    
    Backward.var.new[[T]] = -(lambda[T-1]*tau)^(-2)*((lambda[T-1]*tau)*diag(d)+Sigma_list[[T]])
    
    Backward.Mean.new[T,] = -(lambda[T-1]*tau)^(-1)*mean_matrix[T,]
    
    
    for (t2 in 2:(T-1)){
      
      Forward.var.new[[t2]] = -(lambda[t2]*tau)^(-2)*(solve(Forward.var.new[[t2-1]])+
                                                        (lambda[t2-1]*tau+lambda[t2]*tau)*diag(d)+Sigma_list[[t2]])
      
      Forward.Mean.new[t2,] = -(lambda[t2]*tau)^(-1)*(solve(Forward.var.new[[t2-1]])%*%  Forward.Mean.new[t2-1,]  +
                                                        mean_matrix[t2,])
      
      
      Backward.var.new[[T+1-t2]] = -(lambda[T-t2]*tau)^(-2)*(solve(Backward.var.new[[T+2-t2]])+
                                                               (lambda[T-t2]*tau+lambda[T-t2+1]*tau)*diag(d)+Sigma_list[[T+1-t2]])
      
      Backward.Mean.new[T+1-t2,] = -(lambda[T-t2]*tau)^(-1)*(solve(Backward.var.new[[T+2-t2]])%*% Backward.Mean.new[T+2-t2,] +
                                                               mean_matrix[T+1-t2,])
    }
    
    
    
    X_sigma[[1]] =(solve( Sigma_list[[1]]+solve(Backward.var.new[[2]])+lambda[1]*tau*diag(d)+solve(sigma_X1^2*diag(d))))
    X_sigma[[1]] = (X_sigma[[1]]+t(X_sigma[[1]]))/2
    X_Mean[1,] =  X_sigma[[1]] %*% (mean_matrix[1,]+solve(Backward.var.new[[2]])%*% Backward.Mean.new[2,])
    
    
    
    
    for (t in 2:(T-1)){
      X_sigma[[t]] = (solve(Sigma_list[[t]]+lambda[t-1]*tau*diag(d)+lambda[t]*tau*diag(d)
                           +solve(Forward.var.new[[t-1]])+solve(Backward.var.new[[t+1]])))
      X_sigma[[t]] = (X_sigma[[t]]+t(X_sigma[[t]]))/2
      X_Mean[t,] = X_sigma[[t]]%*%(mean_matrix[t,]+
                                     solve(Forward.var.new[[t-1]])%*%  Forward.Mean.new[t-1,]+
                                     solve(Backward.var.new[[t+1]])%*% Backward.Mean.new[t+1,])
    }
    
    X_sigma[[T]] = solve(Sigma_list[[T]]+lambda[T-1]*tau*diag(d)+solve(Forward.var.new[[T-1]]))
    X_sigma[[T]] = (X_sigma[[T]]+t(X_sigma[[T]]))/2
    X_Mean[T,] =  X_sigma[[T]] %*% (mean_matrix[T,]+solve(Forward.var.new[[T-1]])%*%  Forward.Mean.new[T-1,])
    #   
    
    
    
    Cross_cov[[1]][1:d,1:d] = sigma_X1^(-2)*diag(d)+(lambda[1]*tau)*diag(d)+Sigma_list[[1]]
    Cross_cov[[1]][1:d,(d+1):(2*d)] = -(lambda[1]*tau)*diag(d)
    Cross_cov[[1]][(d+1):(2*d),1:d] = -(lambda[1]*tau)*diag(d)
    Cross_cov[[1]][(d+1):(2*d),(d+1):(2*d)] = (lambda[1]*tau)*diag(d)+solve(Backward.var.new[[3]])+Sigma_list[[2]]+diag(d)+(lambda[2]*tau)
    Cross_cov[[1]] = solve(Cross_cov[[1]])
    
    for ( t in 2:(T-2)){
      Cross_cov[[t]][1:d,1:d] = (lambda[t-1]*tau)*diag(d)+ solve(Forward.var.new[[t-1]]) +Sigma_list[[t]]+(lambda[t]*tau)*diag(d)
      Cross_cov[[t]][1:d,(d+1):(2*d)] = -(lambda[t]*tau)*diag(d)
      Cross_cov[[t]][(d+1):(2*d),1:d] = -(lambda[t]*tau)*diag(d)
      Cross_cov[[t]][(d+1):(2*d),(d+1):(2*d)] = (lambda[t]*tau)*diag(d)+solve(Backward.var.new[[t+2]])+Sigma_list[[t+1]]+(lambda[t+1]*tau)*diag(d)
      Cross_cov[[t]] = solve(Cross_cov[[t]])
    }
    
    Cross_cov[[T-1]][1:d,1:d] = (lambda[T-2]*tau)*diag(d)+ solve(Forward.var.new[[T-2]]) +Sigma_list[[T-1]]+(lambda[T-1]*tau)*diag(d)
    Cross_cov[[T-1]][1:d,(d+1):(2*d)] = -(lambda[T-1]*tau)*diag(d)
    Cross_cov[[T-1]][(d+1):(2*d),1:d] = -(lambda[T-1]*tau)*diag(d)
    Cross_cov[[T-1]][(d+1):(2*d),(d+1):(2*d)] = (lambda[T-1]*tau)*diag(d)+Sigma_list[[T]]
    Cross_cov[[T-1]] = solve(Cross_cov[[T-1]])
    
    
    
    for (t in 1:(T-1)){
      v[t] = 1/(1+lambda[t])
    }
    
    
    b = 0
    b0 = 0
    for (t in 1:(T-1)){
      EX[t] = sum(X_Mean[t+1,]^2)+sum(X_Mean[t,]^2)-2*t(X_Mean[t,]) %*% X_Mean[t+1,]+ 
        sum(diag(x= as.matrix(Cross_cov[[t]][1:d,1:d]+ Cross_cov[[t]][(d+1):(2*d),(d+1):(2*d)] -
                                2*Cross_cov[[t]][1:d,(d+1):(2*d)])))
      
      b = b+ EX[t]*lambda[t]
    }
    b0 = b0 + sum(X_Mean[1,]^2)+sum(diag(x=X_sigma[[1]]))
    
    
    for( t in 1:(T-1)){
      lambda[t] =(d+1)/2/(v[t] + EX[t]*tau/2)
      }
    
    if (global_prior =='Gamma'){
      
      gig_deg = -(T-1)*d/2+1/4
      
      tau = exp(besselK.nuAsym(x=sqrt(b), nu=-gig_deg-1,k.max=5,log=T)-besselK.nuAsym(x=sqrt(b), nu=-gig_deg,k.max=5,log=T))/sqrt(b)-2*gig_deg/b
      
    } else if (global_prior =='Cauthy'){
      v0 = 1/(1+tau)
      
     # tau = (d*(T-1)/2+1/2)/(v0+b/2)
      tau = 1000
    }
    
    
    
     sg_inve_sq = (d+1)/(b0+1)
    # 
     sigma_X1 = sqrt(1/sg_inve_sq)
    
    # lambda_0 = (d/2+1/2)/(1/(1+lambda_0)+b0*lambda_0*tau_0/2)
    # 
    # sigma_X1 = sqrt(1/(tau_0*lambda_0))
    
    # cat(tau,'\n')
    
    
    err[k]=sqrt(sum((X_Mean-Y)^2/T))
    
    
    
    if(abs(err[k]-err[k-1])< gap || k>max_iter || is.nan(tau)){
      ind = 1
    }     else{

      
      X_Mean_old = X_Mean
      
      X_Sigma_old = X_sigma

      k= k+1
    }

  }, silent = TRUE)

if (inherits(res, "try-error")) {
  break
}
}

  X_list= list(Mean=X_Mean_old, Sigma = X_Sigma_old , lambda= 1/sqrt(lambda), v=1/sqrt(v),tau=1/sqrt(tau),err= err[k-1],iter=k, sigma_X1 = sigma_X1)
  return(X_list)
}