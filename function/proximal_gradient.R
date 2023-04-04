# code for proximal gradient comes from  https://rpubs.com/leexiner/bigdata-exercise06


# Negative log likelihood function
nll <- function(X, Y, beta) {
  A <- Y - X %*% beta
  loglike <- (0.5/nrow(X)) * crossprod(A)
  return(loglike)
}
# Gradient of negative log likelihood function
gradient <- function(X, Y, beta) {
  A <- Y - X %*% beta
  grad <- -(1/nrow(X)) * crossprod(X, A)
  return(grad)
}
# Proximal operator
prox.l1 <- function(u, lambda) {
  uhat <- abs(u) - lambda
  ubind <- cbind(rep(0, length(u)), uhat)
  prox <- sign(u) * apply(ubind, 1, max) # 1: row, 2: column
  return(prox)
}
# Proximal gradient descent
proxGD <- function(X,Y,lambda=0.01,gamma=0.01,beta0=NA,iter=1e3,conv=1e-5){
  # Set beta0 equal to a series of zeros
  if (is.na(beta0)) { beta0 <- rnorm(ncol(X)) }
  
  # Initialize coefficients matrix 
  beta <- matrix(rep(NA, ncol(X)*(iter+1)), nrow=ncol(X))
  beta[, 1] <- beta0
  
  # Initialize objective funtion
  obj <- rep(NA, iter+1)
  
  obj[1] <- nll(X, Y, beta0) + lambda * sum(abs(beta0))
  
  # Initialize convergence message in case convergence not reached
  message <- "Convergence not reached..."
  
  for (t in 1:iter) {
    # Update u, beta and obj
    u <- beta[,t] - gamma * gradient(X, Y, beta[,t])
    beta[,t+1] <- prox.l1(u, gamma * lambda)
    obj[t+1] <- nll(X, Y, beta[,t+1]) + lambda * sum(abs(beta[,t+1]))
    
    # Check convergence
    delta <- abs(obj[t+1]-obj[t]) /(abs(obj[t])+conv)
    
    if (delta < conv) {
      # Remove excess betas and obj
      beta <- beta[, -((t+2):ncol(beta))]
      obj <- obj[-((t+2):length(obj))]
      
      # Update convergence message
      message <- sprintf("Convergence reached after %i iterations", (t+1))
      break
    }
  }
  
  result <- list("beta.hat"=beta[,ncol(beta)],"beta"=beta,"objective"=obj,"conv"=message)
  return(result)
}
# Accelerated proximal gradient method
proxACCE <- function(X,Y,lambda=0.01,gamma=0.01,beta0=NA,iter=1e3,conv=1e-5){
  # Set beta0 equal to a series of zeros
  if (is.na(beta0)) { beta0 <- rnorm(ncol(X)) }
  
  # Create s vector
  s <- rep(NA, iter+1)
  s[1] <- 1
  
  for (j in 2:length(s)) {
    s[j] <- (1+sqrt(1+4*s[j-1]^2))/2
  }
  
  # Initialize z matrix
  z <- matrix(0, nrow=ncol(X), ncol=iter+1)
  
  # Initialize coefficients matrix 
  beta <- matrix(rep(NA, ncol(X)*(iter+1)), nrow=ncol(X))
  beta[, 1] <- beta0
  
  # Initialize objective funtion
  obj <- rep(NA, iter+1)
  obj[1] <- nll(X, Y, beta0) + lambda * sum(abs(beta0))
  
  # Initialize convergence message in case convergence not reached
  message <- "Convergence not reached..."
  
  for (t in 1:iter) {
    # Update u, beta, z and obj
    u <- z[,t] - gamma * gradient(X, Y, z[,t])
    beta[,t+1] <- prox.l1(u, gamma * lambda)
    z[,t+1] <- beta[,t+1] + (s[t] - 1)/s[t+1] * (beta[,t+1] - beta[,t])
    obj[t+1] <- nll(X, Y, beta[,t+1]) + lambda * sum(abs(beta[,t+1]))
    
    # Check convergence 
    delta <- abs(obj[t+1]-obj[t]) /(abs(obj[t])+conv)
  
    
    if (delta < conv) {
      # Remove excess betas and nll
      beta <- beta[, -((t+2):ncol(beta))]
      obj <- obj[-((t+2):length(obj))]
      
      # Update convergence message
      message <- sprintf("Convergence reached after %i iterations", (t+1))
      break
    }
  }
  
  result <- list("beta.hat"=beta[,ncol(beta)],"beta"=beta,"objective"=obj,"conv"=message)
  return(result)
}