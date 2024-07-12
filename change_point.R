 rm(list=ls()) 
 set.seed(2024)
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
 
 
# Initialize CSV file with column names
 write.table(data.frame(
   method = character(),
   iterations = integer(),
   rho = numeric(),
   error = numeric()
 ), file = "change_point_simu.csv", append = FALSE, sep = ",", col.names = TRUE, row.names = FALSE)

 
 iter_max = 50
 

 for (rho in c(1,0.75,0.5,0.25)){  # sparsity parameter
for (iter in 1:iter_max){
  results <- data.frame(
    method = character(),
    iterations = integer(),
    rho = numeric(),
    error = numeric()
  )
  
 
 
 p = 15 # number of nodes
 
 block_num = 3 # number of groups for SBM
 n = 50 # sample size for each segment
 # connectivity matrix for the first and the third segments
 conn1_mat =  rho*matrix(c(0,1,0,1,0,0,0,0,1), nrow = 3)
 # connectivity matrix for the second segment
 conn2_mat =  rho*matrix(c(1,0,0,0,0,1,0,1,0), nrow = 3)
 can_vec = sample(1:p, replace = FALSE) # randomly assign nodes into groups
 sbm1 = simu.SBM(conn1_mat, can_vec, n, symm = TRUE, self = TRUE)
 sbm2 = simu.SBM(conn2_mat, can_vec, n, symm = TRUE, self = TRUE)
 data_mat = cbind(sbm1$obs_mat, sbm2$obs_mat)
 ### detect change points
 M = 10 # number of random intervals for WBS
 d = 5 # parameter for scaled PCA algorithm
 delta = 5
 intervals = WBS.intervals(M = M, lower = 1, upper = ncol(data_mat))
 WBS_result = WBS.nonpar.RDPG(data_mat, lowerdiag = TRUE,  d,
                              Alpha = intervals$Alpha, Beta = intervals$Beta, delta)
 # cpt_hat = tuneBSnonparRDPG(WBS_result, data_mat, lowerdiag = TRUE, d)
 # cpt_hat 

 
 T = 2*n
 Y = vector("list", T)
 
 for (t in 1:T){
   Y[[t]] =  lowertri2mat(data_mat[,t], p, diag = FALSE)
 }
 
 MF_list = MP_binary_weighted_adaptive (Y, gap =1e-3, max_iter=200,
                                        global_prior='Cauthy',first_index =TRUE)
 
 FFS_result = FFS.RDPG(xhat =MF_list$Mean_X,betahat=MF_list$mean_beta, obs_num=T,
                              Alpha = intervals$Alpha, Beta = intervals$Beta, delta=5)
 
 
 
 error = abs(WBS_result$S[which.max(WBS_result$Dval)]-50)
 
 cat('WBS',error,'\n')
 # Append results
 results <- rbind(results, data.frame(
   method = "WBS_RDP",
   iterations = iter,
   rho = rho,
   error = error
 ))
 
 error = abs(FFS_result$S[which.max(FFS_result$Dval)]-50)
 
 cat('FFS',error,'\n')
 # Append results
 results <- rbind(results, data.frame(
   method = "WBS_FFS",
   iterations = iter,
   rho = rho,
   error = error
 ))
 write.table(results, file = "change_point_simu.csv", append = TRUE, sep = ",", col.names = FALSE, row.names = FALSE)
 
}
   
   
}
 
 
 
# output = cbind(FFS_err,WBS_err)
# 
# write.table(output, file = "change_point.csv", append = TRUE, sep = ",", col.names = FALSE, row.names = FALSE)
# 
# 
# 
change_point_df = read.csv('change_point_simu3.csv')

# change_point_df = change_point_df[change_point_df$method!='WBS_RDP',]

change_point_df$n = 'n=15'

change_point_df_2 = read.csv('change_point_simu5.csv')

change_point_df_2$n ='n=30'

change_point_df = rbind(change_point_df,change_point_df_2)

change_point_df$rho = as.character(change_point_df$rho)

change_point_df$method[change_point_df$method=='WBS_FFS']='FFS'

change_point_df$method[change_point_df$method=='WBS_RDP']='RDP'

# 
# change_point_df$error = c(FFS_err[1:200],WBS_err[1:200])
# change_point_df$method = c(rep('FFS',200),rep('WBS',200))
# change_point_df$rho=c(rep(c(rep('1',50),rep('0.75',50),rep('0.5',50),rep('0.25',50)),2))
# 
 change_point_df = as.data.frame(change_point_df)
# 
 library(ggplot2)
 ggplot(change_point_df, aes(x = rho, y = error, fill = method)) +
   geom_boxplot(outlier.shape = NA)+ facet_wrap(vars(factor(n)),scale='free')+ theme(axis.text.x = element_text(colour = 'black', size = 15),
                                                                                                                      axis.title.x = element_text(size = 15,
                                                                                                                                                  hjust = 0.5, vjust = 0.2)) +
   theme(axis.text.y = element_text(colour = 'black', size = 15),
         axis.title.y = element_text(size = 15,
                                     hjust = 0.5, vjust = 0.2)) +
   theme(legend.text=element_text(size=15),legend.title=element_text(size=15),
         legend.key.size = unit(1, 'cm'))+
   theme(strip.text = element_text(size=15))+
   # geom_boxplot(outlier.shape = NA) +
   # coord_cartesian(ylim = c(0,50))+
   # theme(plot.title = element_text(hjust = 0.5))+
   labs(x = expression(Sparsity~factor~rho), y = "Detection errors")