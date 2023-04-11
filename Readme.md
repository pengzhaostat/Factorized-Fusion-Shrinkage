Code for Factorized Fusion Shrinkage for Dynamic Relational Data

### Files for simulation

`motivating_example2.R`: code for reproducing Fig. 1 and Fig. 3

`Gaussian_simu.R`: Example of case 1 in simulation 

`binary_simu.R`: Example of case 2 in simulation 

`tensor_simu.R`:  Example of case 3 in simulation 

`trend_simu_algorithm_compare.R`: code for reproducing Fig. 10

`formal_alliance.R`: code for reproducing formal alliance data analysis, Fig. 7-8 (Up to some orthogonal transformations)

`airline.R`: code for reproducing airline data analysis, Fig. 6

`global_local_prior.R`: code for reproducing Fig. 9

### Files for functions, in the function folder

`MP_gibbs_multi_Sigma_adaptive.R`: core function, implement main VI algorithm for each block

`MP_gibbs_Gaussnetwork_adaptive.R`: function for the Gaussian matrix factorization case

`MP_gibbs_network_adaptive.R`: function for the binary matrix factorization case

`MP_gibbs_network_sym_adaptive.R`: function for the binary network case

`MP_gibbs_tensor_adaptive.R`:  function for the Gaussian tensor decomposition case

`helper.R`: some helper functions

`mix_DN_adaptive_invgamma.R`: comparing method, inverse gamma prior to the transition variances 

`proximal_gradient.R`: comparing method, proximal gradient descent for fused lasso 
