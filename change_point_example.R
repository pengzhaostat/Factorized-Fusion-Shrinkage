rm(list=ls()) 
set.seed(2023)
library(changepoints)
library(MASS)
library(mvtnorm)
library(igraph)
library(DiagrammeR)
library(pROC)
library(Bessel)
library(animation)
file_location = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(paste(file_location,'/function',sep=""))
source('MP_gibbs_multi_Sigma_adaptive.R')
source('MP_gibbs_network_sym_adaptive.R')
source('mix_DN_adaptive_invgamma.R')
source('helper.R')
source('change_point_WB_FFS.R')
setwd(file_location)

p = 20 # number of nodes

rho = 0.5 #signal strength

block_num = 3 # number of groups for SBM
n = 50 # sample size for each segment, true change point t0 = 50
# connectivity matrix for the first and the third segments
conn1_mat =   rho*matrix(c(0,1,0,1,0,0,0,0,1), nrow = 3)
# connectivity matrix for the second segment
conn2_mat =   rho*matrix(c(1,0,0,0,0,1,0,1,0), nrow = 3)
can_vec = sample(1:p, replace = FALSE) # randomly assign nodes into groups
sbm1 = simu.SBM(conn1_mat, can_vec, n, symm = TRUE, self = TRUE)
sbm2 = simu.SBM(conn2_mat, can_vec, n, symm = TRUE, self = TRUE)
data_mat = cbind(sbm1$obs_mat, sbm2$obs_mat)
### detect change points
M = 10 # number of random intervals for WBS
intervals = WBS.intervals(M = M, lower = 1, upper = ncol(data_mat))
d = 5 # parameter for scaled PCA algorithm
delta = 5





T = 2*n
Y = vector("list", T)

for (t in 1:T){
  Y[[t]] =  lowertri2mat(data_mat[,t], p, diag = FALSE)
}

MF_list = MP_binary_weighted_adaptive (Y, gap =1e-3, max_iter=200, d=d,
                                       global_prior='Cauthy',first_index =TRUE)

FFS_result = FFS.RDPG(xhat =MF_list$Mean_X,betahat=MF_list$mean_beta, obs_num=T,
                      Alpha = intervals$Alpha, Beta = intervals$Beta, delta=5)

###RDP using changepoints package
WBS_result = WBS.nonpar.RDPG(data_mat, lowerdiag = TRUE, d,
                             Alpha = intervals$Alpha, Beta = intervals$Beta, delta=5)


WBS_result$S[which.max(WBS_result$Dval)]

FFS_result$S[which.max(FFS_result$Dval)]



