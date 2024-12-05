## Code for Factorized Fusion Shrinkage for Dynamic Relational Data

# Files for simulation
motivating_example2.R: code for reproducing Fig. 1

main_simu_cluster.R:  code for reproducing simulation of Fig.2

main_simu_time.R: code for reproducing simulation of Fig.3

main_simu_change_point.R: code for reproducing simulation of Fig.4

global_local_prior.R: code for reproducing Fig. 5

formal_alliance.R: code for reproducing formal alliance data analysis, Fig. 4 (Up to some orthogonal transformations)

change_point_exmaple.R: Example of simulation for change point detection

cluster_example: Example of simulation for cluster

simu_smooth.R: code for reproducing Fig. 8

Gaussian_simu_dim.R: code for reproducing Fig. 9

Gaussian_simu_hyper.R: code for reproducing Fig. 10

Gaussian_simu.R: Example of case 1 in simulation in the appendix

binary_simu.R: Example of case 2 in simulation in the appendix

tensor_simu.R: Example of case 3 in simulation in the appendix

trend_simu_algorithm_compare.R: code for reproducing Fig. 11

# Files for functions, in the function folder

MP_gibbs_network_sym_adaptive_cpp.cpp: core function, implement main VI algorithm for each block for dynamic networks.

MP_gibbs_multi_Sigma_adaptive.R: R version of core function, which could be slow, implement main VI algorithm for each block

MP_gibbs_Gaussnetwork_adaptive.R: R version of function, which could be slow, for the Gaussian matrix factorization case

MP_gibbs_network_adaptive.R: R version of function for the binary matrix factorization case

MP_gibbs_network_sym_adaptive.R: R version of function for the binary network case

MP_gibbs_tensor_adaptive.R: R version of function for the Gaussian tensor decomposition case

helper.R: some helper functions

mix_DN_adaptive_inv_gamma_new.cpp: comparing method, inverse gamma prior to the transition variances

mix_DN_adaptive_invgamma.R: R version of the comparing method, inverse gamma prior to the transition variances

proximal_gradient.R: comparing method, proximal gradient descent for fused lasso

change_point_WB_FFS: implement the change point detection for dynamic networks after using FFS for denoising

mcmc_DN_adaptive: MCMC algorithm for FFS
