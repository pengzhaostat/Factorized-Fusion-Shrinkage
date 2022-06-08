Code for Factorized Fusion Shrinkage for Dynamic Relational Data

### Files for simulation

`Gaussian_simu.R`: Example of case 1 in simulation 

`binary_simu.R`: Example of case 2 in simulation 

`tensor_simu.R`:  Example of case 3 in simulation 

`motivating_example.R`: code for reproducing Fig. 1

`motivating_example2.R`: code for reproducing Fig. 2

`trend_simu_algorithm_compare.R`: code for reproducing Fig. 2

`formal_alliance.R`: code for reproducing formal alliance data analysis, Fig. 3-5 (Up to some orthogonal transformations)

`animation_alliance.pdf`: animated version of visulization of formal alliance data analysis


### Files for functions, in the function folder

`MP_gibbs_multi_Sigma_adaptive.R`: core function, implement main VI algorithm for each block

`MP_gibbs_Gaussnetwork_adaptive.R`: function for the Gaussian matrix factorization case

`MP_gibbs_network_adaptive.R`: function for the binary matrix factorization case

`MP_gibbs_network_sym_adaptive.R`: function for the binary network case

`MP_gibbs_tensor_adaptive.R`:  function for the Gaussian tensor decomposition case

`helper.R`: some helper functions

`mix_DN_adaptive_invgamma.R`: comparing methood, inverse gamma prior to the transition variances 

`proximal_gradient.R`: comparing methood, proximal gradient descent for fused lasso 
