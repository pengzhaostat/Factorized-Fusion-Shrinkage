FFS.RDPG = function (xhat,betahat, obs_num, lowerdiag = FALSE,  Alpha, Beta, delta) 
{

  n = nrow(xhat[[1]])

  Y_mat = matrix(0, n*(n-1)/2, obs_num)
  for (t in 1:obs_num) {
    phat = drop(1/(1+exp(-betahat-xhat[[t]] %*% t(xhat[[t]]))))
    Y_mat[,t] = phat[lower.tri(phat)]
    # ind = sample(1:n, n, replace = FALSE)
    # for (i in 1:nrow(Y_mat)) {
    #   Y_mat[i, t] = phat[ind[floor(i/n)+1], ind[floor(i/n)+i%%n]]
    # }
  }
  result = WBS.uni.nonpar_FFS(Y_mat, s = 1, e = obs_num, Alpha, 
                          Beta, N = rep(nrow(Y_mat),obs_num), delta)
  return(result)
}

CUSUM.KS = function(Y, s, e, t, N, vector = FALSE){
  n_st = log(sum(N[s:t]))
  n_se =  log(sum(N[s:e]))
  n_te =  log(sum(N[(t+1):e]))
  aux = as.vector(Y[,s:t])
  aux = aux[which(is.na(aux)==FALSE)]
  temp = ecdf(aux)
  vec_y = as.vector(Y[,s:e])
  vec_y = vec_y[which(is.na(vec_y)==FALSE)]
  Fhat_st = temp(vec_y)# temp(grid)
  aux = as.vector(Y[,(t+1):e])
  aux = aux[which(is.na(aux)==FALSE)]
  temp = ecdf(aux)
  Fhat_te = temp(vec_y)# temp(grid)
  if(vector == TRUE){
    result = sqrt(exp(n_st + n_te - n_se)) * abs(Fhat_te - Fhat_st)
  }else{
    result = sqrt(exp(n_st + n_te - n_se)) * max(abs(Fhat_te - Fhat_st)) 
  }
  return(result)
}

WBS.uni.nonpar_FFS = function (Y, s, e, Alpha, Beta, N, delta, level = 0) 
{
  Alpha_new = pmax(Alpha, s)
  Beta_new = pmin(Beta, e)
  idx = which(Beta_new - Alpha_new > 2 * delta)
  Alpha_new = Alpha_new[idx]
  Beta_new = Beta_new[idx]
  M = length(Alpha_new)
  S = NULL
  Dval = NULL
  Level = NULL
  Parent = NULL
  if (M == 0) {
    return(list(S = S, Dval = Dval, Level = Level, Parent = Parent))
  }  else {
    level = level + 1
    parent = matrix(c(s, e), nrow = 2)
    a = rep(0, M)
    b = rep(0, M)
    for (m in 1:M) {
      temp = rep(0, Beta_new[m] - Alpha_new[m] - 2 * delta + 
                   1)
      for (t in (Alpha_new[m] + delta):(Beta_new[m] - delta)) {
        temp[t - (Alpha_new[m] + delta) + 1] = CUSUM.KS(Y,Alpha_new[m], Beta_new[m], t, N)
      }
      best_value = max(temp)
      best_t = which.max(temp) + Alpha_new[m] + delta -  1
      a[m] = best_value
      b[m] = best_t
    }
    m_star = which.max(a)
  }
  temp1 = WBS.uni.nonpar_FFS(Y, s, b[m_star] - 1, Alpha, Beta, 
                         N, delta, level)
  temp2 = WBS.uni.nonpar_FFS(Y, b[m_star], e, Alpha, Beta, N, delta, 
                         level)
  S = c(temp1$S, b[m_star], temp2$S)
  Dval = c(temp1$Dval, a[m_star], temp2$Dval)
  Level = c(temp1$Level, level, temp2$Level)
  Parent = cbind(temp1$Parent, parent, temp2$Parent)
  result = list(S = S, Dval = Dval, Level = Level, Parent = Parent)
  class(result) = "BS"
  return(result)
}