#define ARMA_64BIT_WORD 1

#include <R.h>
#include <Rmath.h>

// [[Rcpp::plugins(cpp11)]]
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <cmath>

#include <vector>
using std::vector;

#include <iostream>
using std::cout;
using std::endl;

using namespace Rcpp;
using namespace arma;

/* supporting function */
// [[Rcpp::export]]
double distance_squared_inner_prod(arma::vec mu1, arma::vec mu2, arma::mat Sigma1, arma::mat Sigma2){
  // mu1 and mu2 are d x 1, Sigma1 and Sigma2 are d x d
  
  double S = 0.0;
  double mu12 = dot(mu1, mu2);
  S = mu12 * mu12 + sum((Sigma1 * Sigma2).eval().diag()) + (mu1.t() * Sigma2 * mu1).eval()(0,0) + (mu2.t() * Sigma1 * mu2).eval()(0,0);
  return S;
}

Rcpp::List MP_gibbs_mult_Sigma(const mat& Y,  
                               const cube& Sigma_list ,double tau =1, int max_iter=200, 
                               double gap = 1e-3) {
  
  int T = Y.n_rows;
  int d = Y.n_cols;
  int K = 2000;
  
  mat X_Mean = arma::randn<mat>(T, d);
  mat Forward_Mean_new = arma::zeros<mat>(T, d);
  mat Backward_Mean_new = arma::zeros<mat>(T, d);
  cube X_sigma(d, d, T);
  cube Forward_var_new(d, d, T);
  cube Backward_var_new(d, d, T);
  cube Cross_cov(2 * d, 2 * d, T);
  
  mat mean_matrix = Y;
  
  for (int t = 0; t < T; ++t) {
    X_sigma.slice(t) = eye<mat>(d, d);
    Forward_var_new.slice(t) = eye<mat>(d, d);
    Backward_var_new.slice(t) = eye<mat>(d, d);
    Cross_cov.slice(t) = eye<mat>(2 * d, 2 * d);
    mean_matrix.row(t) = (Sigma_list.slice(t) * Y.row(t).t()).t();
  }
  
  vec lambda(T - 1, fill::ones);
  vec lambda_new(T - 1, fill::ones);
  vec v(T - 1, fill::ones);
  
  vec EX(T - 1, fill::zeros);
  vec err(K, fill::zeros);
  
  err(0) = sqrt(sum(sum(pow(X_Mean - Y, 2)) / T));
  
  vec err_e;
  

  double v0 = 1;
  int k = 1;
  int ind = 0;
  
  double sigma_X1 = 0.5;
  
  mat X_Mean_old = X_Mean;
  cube X_Sigma_old = X_sigma;
  
  
  
  while (k < K && ind == 0) {
    
    Forward_var_new.slice(0) = -(pow(lambda(0) * tau, -2)) * (pow(sigma_X1, -2) * eye<mat>(d, d) + lambda(0) * tau * eye<mat>(d, d) + Sigma_list.slice(0));
    Forward_Mean_new.row(0) = -(pow(lambda(0) * tau, -1)) * mean_matrix.row(0);
    
    
    Backward_var_new.slice(T - 1) = -(pow(lambda(T - 2) * tau, -2)) * ((lambda(T - 2) * tau) * eye<mat>(d, d) + Sigma_list.slice(T - 1));
    Backward_Mean_new.row(T - 1) = -(pow(lambda(T - 2) * tau, -1)) * mean_matrix.row(T - 1);
    
    for (int t2 = 1; t2 < (T - 1); ++t2) {
      Forward_var_new.slice(t2) = -(pow(lambda(t2) * tau, -2)) * (inv(Forward_var_new.slice(t2 - 1)) +
        (lambda(t2 - 1) * tau + lambda(t2) * tau) * eye<mat>(d, d) + Sigma_list.slice(t2));
      Forward_Mean_new.row(t2) = -(pow(lambda(t2) * tau, -1)) * ( (inv(Forward_var_new.slice(t2 - 1)) * Forward_Mean_new.row(t2 - 1).t()).t() +
        mean_matrix.row(t2));
      
      Backward_var_new.slice(T - 1 - t2) = -(pow(lambda(T - 2 - t2) * tau, -2)) * (inv(Backward_var_new.slice(T - t2)) +
        (lambda(T - 2 - t2) * tau + lambda(T - 1 - t2) * tau) * eye<mat>(d, d) + Sigma_list.slice(T - 1 - t2));
      Backward_Mean_new.row(T - 1 - t2) = -(pow(lambda(T - 2 - t2) * tau, -1)) * ((inv(Backward_var_new.slice(T - t2)) * Backward_Mean_new.row(T - t2).t()).t() +
        mean_matrix.row(T - 1 - t2));
      
    }
    
    X_sigma.slice(0) = inv(Sigma_list.slice(0) + inv(Backward_var_new.slice(1)) + lambda(0) * tau * eye<mat>(d, d) + inv(pow(sigma_X1, 2) * eye<mat>(d, d)));
    
    X_sigma.slice(0) = (X_sigma.slice(0)+X_sigma.slice(0).t())/2;
    
    
    X_Mean.row(0) = (X_sigma.slice(0) * (mean_matrix.row(0).t() + inv(Backward_var_new.slice(1)) * Backward_Mean_new.row(1).t())).t();
    
    
    for (int t = 1; t < (T - 1); ++t) {
      X_sigma.slice(t) = inv(Sigma_list.slice(t) + lambda(t - 1) * tau * eye<mat>(d, d) + lambda(t) * tau * eye<mat>(d, d) +
        inv(Forward_var_new.slice(t - 1)) + inv(Backward_var_new.slice(t + 1)));
      X_sigma.slice(t) = (X_sigma.slice(t)+X_sigma.slice(t).t())/2;
      
      X_Mean.row(t) = (X_sigma.slice(t) * (mean_matrix.row(t).t() +
        (inv(Forward_var_new.slice(t - 1)) * Forward_Mean_new.row(t - 1).t()) +
        (inv(Backward_var_new.slice(t + 1)) * Backward_Mean_new.row(t + 1).t()))).t();
    }
    
    X_sigma.slice(T - 1) = inv(Sigma_list.slice(T - 1) + lambda(T - 2) * tau * eye<mat>(d, d) + inv(Forward_var_new.slice(T - 2)));
    
    X_sigma.slice(T - 1) = (X_sigma.slice(T - 1)+X_sigma.slice(T - 1).t())/2;
    
    X_Mean.row(T - 1) = (X_sigma.slice(T - 1) * (mean_matrix.row(T - 1).t() +
      (inv(Forward_var_new.slice(T - 2)) * Forward_Mean_new.row(T - 2).t()))).t();
    
    Cross_cov.slice(0).submat(0, 0, d - 1, d - 1) = pow(sigma_X1, -2) * eye<mat>(d, d) + lambda(0) * tau * eye<mat>(d, d) + Sigma_list.slice(0);
    Cross_cov.slice(0).submat(0, d, d - 1, 2 * d - 1) = -lambda(0) * tau * eye<mat>(d, d);
    Cross_cov.slice(0).submat(d, 0, 2 * d - 1, d - 1) = -lambda(0) * tau * eye<mat>(d, d);
    Cross_cov.slice(0).submat(d, d, 2 * d - 1, 2 * d - 1) = lambda(0) * tau * eye<mat>(d, d) + inv(Backward_var_new.slice(2)) +
      Sigma_list.slice(1) +  lambda(1) * tau * eye<mat>(d, d);
    Cross_cov.slice(0) = inv(Cross_cov.slice(0));
    
    for (int t = 1; t < (T - 2); ++t) {
      Cross_cov.slice(t).submat(0, 0, d - 1, d - 1) = lambda(t - 1) * tau * eye<mat>(d, d) + inv(Forward_var_new.slice(t - 1)) +
        Sigma_list.slice(t) + lambda(t) * tau * eye<mat>(d, d);
      Cross_cov.slice(t).submat(0, d, d - 1, 2 * d - 1) = -lambda(t) * tau * eye<mat>(d, d);
      Cross_cov.slice(t).submat(d, 0, 2 * d - 1, d - 1) = -lambda(t) * tau * eye<mat>(d, d);
      Cross_cov.slice(t).submat(d, d, 2 * d - 1, 2 * d - 1) = lambda(t) * tau * eye<mat>(d, d) +
        inv(Backward_var_new.slice(t + 2)) + Sigma_list.slice(t + 1) +
        lambda(t + 1) * tau * eye<mat>(d, d);
      Cross_cov.slice(t) = inv(Cross_cov.slice(t));
    }
    
    Cross_cov.slice(T - 2).submat(0, 0, d - 1, d - 1) = lambda(T - 3) * tau * eye<mat>(d, d) + inv(Forward_var_new.slice(T - 3)) +
      Sigma_list.slice(T - 2) + lambda(T - 2) * tau * eye<mat>(d, d);
    Cross_cov.slice(T - 2).submat(0, d, d - 1, 2 * d - 1) = -lambda(T - 2) * tau * eye<mat>(d, d);
    Cross_cov.slice(T - 2).submat(d, 0, 2 * d - 1, d - 1) = -lambda(T - 2) * tau * eye<mat>(d, d);
    Cross_cov.slice(T - 2).submat(d, d, 2 * d - 1, 2 * d - 1) = lambda(T - 2) * tau * eye<mat>(d, d) + Sigma_list.slice(T - 1);
    Cross_cov.slice(T - 2) = inv(Cross_cov.slice(T - 2));
    
    for (int t = 0; t < (T - 1); ++t) {
      v(t) = 1 / (1 + lambda(t));
    }
    
    double b = 0;
    double b0 = 0;
    for (int t = 0; t < (T - 1); ++t) {
      EX[t] = dot(X_Mean.row(t + 1), X_Mean.row(t + 1)) + dot(X_Mean.row(t), X_Mean.row(t)) - 
        2 * dot(X_Mean.row(t), X_Mean.row(t + 1)) +
        sum(diagvec(Cross_cov.slice(t).submat(0, 0, d - 1, d - 1) +
        Cross_cov.slice(t).submat(d, d, 2 * d - 1, 2 * d - 1) -
        2 * Cross_cov.slice(t).submat(0, d, d - 1, 2 * d - 1)));
      b += EX[t] * lambda[t];
    }
    b0 += dot(X_Mean.row(0),X_Mean.row(0)) + sum((X_sigma.slice(0)).diag());
    
    for (int t = 0; t < (T - 1); ++t) {
      lambda[t] = (d + 1) / (2 * (v[t] + EX[t] * tau / 2));
    }
    
    v0 = 1/(1+tau);
    
    tau = (d*(T-1)/2+1/2)/(v0+b/2);
    
    
    sigma_X1 = sqrt(1 / ((d + 1) / (b0 + 1)));
    
    err(k) = sqrt(sum(sum(pow(X_Mean - Y, 2)) / T));
    
    if (std::abs(err(k) - err(k - 1)) < gap || k > max_iter || std::isnan(tau)) {
      ind = 1;
    } else {
      X_Mean_old = X_Mean;
      X_Sigma_old = X_sigma;
      ++k;
    }
  }
  
  return  Rcpp::List::create(
    Rcpp::Named("Mean") = X_Mean_old,
    Rcpp::Named("Sigma") = X_Sigma_old,
    Rcpp::Named("lambda") = 1 / sqrt(lambda),
    Rcpp::Named("v") = 1 / sqrt(v),
    Rcpp::Named("tau") = 1 / sqrt(tau),
    Rcpp::Named("err") = err,
    Rcpp::Named("iter") = k,
    Rcpp::Named("sigma_X1") = sigma_X1);
}
// Define the main function MP_binary_weighted_adaptive
// [[Rcpp::export]]
Rcpp::List MP_binary_weighted_adaptive(const Cube<double>& Y,
                                 double tau = 0.01,
                                 double gap = 1e-4,
                                 int max_iter = 2000,
                                 double alpha = 0.95,
                                 int d = 2,
                                 double mean_beta_prior = 0,
                                 double inter_gap = 1e-4,
                                 double sigma_beta_prior = 3) {
  
  // Initialization
  tau = 1 / pow(tau, 2);
  int T = Y.n_slices;
  int n = Y.n_rows;
  
  arma::cube Mean_X(n, d, T, fill::zeros);
  Mean_X.randn();
//  Mean_X *= 0.1;
  
  std::vector<cube> Sigma_X(n);
  for (int i = 0; i < n; ++i) {
    Sigma_X[i] = cube(d, d, T, fill::zeros);
  }
  
  cube Mean_X_new = Mean_X;
  std::vector<cube> Sigma_X_new(n);
  for (int i = 0; i < n; ++i) {
    Sigma_X_new[i] = cube(d, d, T, fill::zeros);
  }
  
  mat lambda_X(n, T - 1, fill::ones);
  mat v_X(n, T - 1, fill::ones);
  
  double tau_X;
  tau_X =tau_X + tau;
  vec EX(T - 1, fill::zeros);
  
  vec M_b(d, fill::zeros);
  mat V_b(d, d, fill::zeros);
  
  arma::cube Xi(n, n, T, fill::zeros);
  
  int K = 1000;
  vec err(K, fill::zeros);
  err(0) = 0;
  vec norm_stop(K, fill::zeros);
  int k = 2;
  int ind = 0;
  vec pred_mean(T * n * n, fill::zeros);
  vec pred_mean_old(T * n * n, fill::zeros);
  
  arma::vec vec1( (n-1)*n*T, arma::fill::zeros);
  
  mat  M_X(T,d, fill::zeros);
  cube V_X(n,d,T, fill::zeros);
  
  double sigma_beta = 0;
  
  double mean_beta = 0;
  
  while (k < K && ind == 0) {
    double V_beta_cumulative = 0;
    double M_beta_cumulative = 0;
    
    std::vector<cube> V_X_cumulative(n);
    arma::cube M_X_cumulative(n, d, T, fill::zeros);
    
    for (int i = 0; i < n; ++i) {
      V_X_cumulative[i] = cube(d, d, T, fill::zeros);
    }
    
    // Update of beta
    for (int t = 0; t < T; ++t) {
      for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
          vec M_i = Mean_X.slice(t).row(i).t();
          vec M_j = Mean_X.slice(t).row(j).t();
          mat V_i = Sigma_X[i].slice(t);
          mat V_j = Sigma_X[j].slice(t);
          double c = sqrt(distance_squared_inner_prod(M_i, M_j, V_i, V_j) + pow(mean_beta, 2) + pow(sigma_beta, 2));
          Xi(i, j, t) = 1 / (2 * c) * (1 - exp(-c)) / (exp(-c) + 1) / (-2);
          if (j != i) {
            V_beta_cumulative -= 2 * Xi(i, j, t) * alpha;
            M_beta_cumulative += (Y(i, j, t) - 0.5 + 2 * Xi(i, j, t) * dot(M_i, M_j)) * alpha;
          }
        }
      }
      
    }
    

    // Update of sigma_beta
    double sigma_beta_new = 1/sqrt(pow(sigma_beta_prior, -2) + V_beta_cumulative);
    
    // Update of mean_beta
    double mean_beta_new = sigma_beta_new * sigma_beta_new * (pow(sigma_beta_prior, -2) * mean_beta_prior + M_beta_cumulative);
    
    
    // Update of Sigma_X and Mean_X
    arma::uvec subset_indices = arma::randperm(n);
    for (unsigned int i0 = 0; i0 < subset_indices.n_elem; i0++) {
      unsigned int i = subset_indices(i0);
      for (int t = 0; t < T; ++t) {
        for (unsigned int j0 = 0; j0 < subset_indices.n_elem; j0++) {
          unsigned int j = subset_indices(j0);
          if (j0 > i0) {
            vec M_j = Mean_X.slice(t).row(j).t();
            mat V_j = Sigma_X[j].slice(t);
            V_X_cumulative[i].slice(t) -= 2 * Xi(i, j, t) * (M_j * M_j.t() + V_j) * alpha;
            M_X_cumulative.slice(t).row(i) += (Y(i, j, t) - 0.5 + 2 * Xi(i, j, t) * mean_beta_new) * M_j.t() * alpha;
          } else if (j0 < i0) {
            vec M_j = Mean_X_new.slice(t).row(j).t();
            mat V_j = Sigma_X_new[j].slice(t);
            V_X_cumulative[i].slice(t) -= 2 * Xi(i, j, t) * (M_j * M_j.t() + V_j) * alpha;
            M_X_cumulative.slice(t).row(i) += (Y(i, j, t) - 0.5 + 2 * Xi(i, j, t) * mean_beta_new) * M_j.t() * alpha;
          }
        }

        
        mat V_b = V_X_cumulative[i].slice(t);
       
        vec M_b = M_X_cumulative.slice(t).row(i).t();
  
    
        M_X_cumulative.slice(t).row(i) = (inv(V_b) * M_b).t();

      }
      mat M_X_i =  M_X_cumulative.row(i);
      cube V_X_i =  V_X_cumulative[i];
 
      // Update Sigma_X and Mean_X using MP_gibbs_mult_Sigma function
      Rcpp::List MF_gibbs =  MP_gibbs_mult_Sigma (M_X_i.t(),  V_X_i, tau=tau_X, 100, inter_gap);
      
      // Assign updated Mean_X and Sigma_X_new values
      M_X = Rcpp::as<arma::mat>(MF_gibbs["Mean"]);
      V_X = Rcpp::as<arma::cube>(MF_gibbs["Sigma"]);
 //     tau_X = Rcpp::as<double>(MF_gibbs["tau"]);
      Mean_X_new.row(i) = M_X.t();
      Sigma_X_new[i] = V_X;
//      lambda_X(i) = Rcpp::as<arma::vec>(MF_gibbs["lambda"]);


    }
    
    
    
    vec auc_mean = vec1;
    vec rmse_iter = vec1;
    int r = 0;
  for(int t=0; t<T; t++){
    for( int i=0; i<n; i++){
      for( int j=0; j<n; j++){
        if(j < i){
          auc_mean(r) = 1/(1+exp( -mean_beta_new -arma::dot(Mean_X_new.slice(t).row(i), Mean_X_new.slice(t).row(i)) ));
          rmse_iter(r) = 1/(1+exp( -mean_beta -arma::dot( Mean_X.slice(t).row(i), Mean_X.slice(t).row(i) ) ));
          r = r + 1;
        }
      }
    }
    }
    
    double temp_s = arma::dot(auc_mean-rmse_iter, auc_mean-rmse_iter);
    err(k-1) = sqrt(temp_s/(n-1)/n/T*2);
    
    if( (abs(err(k-1))<gap) || (k>max_iter) ){ ind = 1; }
    else{
      printf("Iteration %d: %f\n", k, err(k-1));
      Mean_X = Mean_X_new;

      Sigma_X = Sigma_X_new;
      mean_beta = mean_beta_new;
      sigma_beta = sigma_beta_new;
      k = k + 1;
    }
  }
  
  return Rcpp::List::create(Named("Mean_X") = Mean_X,
                      Named("Sigma_X") = Sigma_X,
                      Named("iter") = k - 1,
                      Named("tau") = tau_X,
                      Named("lambda_X") = lambda_X,
                      Named("mean_beta") = mean_beta,
                      Named("sigma_beta") = sigma_beta,
                      Named("Xi") = Xi,
                    Named("err") = err);
}
