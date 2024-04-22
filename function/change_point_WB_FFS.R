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
  result = WBS.uni.nonpar(Y_mat, s = 1, e = obs_num, Alpha, 
                          Beta, N = rep(nrow(Y_mat), obs_num), delta)
  return(result)
}
