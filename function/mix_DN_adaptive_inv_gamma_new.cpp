// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace arma;
using namespace Rcpp;

// Function to compute distance squared inner product
double distance_squared_inner_prod(const arma::vec& mu1, const arma::vec& mu2, const arma::mat& Sigma1, const arma::mat& Sigma2) {
  double mu12 = dot(mu1, mu2);
  double S = mu12 * mu12 + trace(Sigma1 * Sigma2) + as_scalar(mu1.t() * Sigma2 * mu1) + as_scalar(mu2.t() * Sigma1 * mu2);
  return S;
}

// [[Rcpp::export]]
Rcpp::List mix_DN_adaptive_invgamma(const std::vector<arma::mat>& Y, // Y is a vector of matrices
                                    int d = 2,
                                    double mean_beta_prior = 0,
                                    double sigma_beta_prior = 3,
                                    double gap = 0.01,
                                    int min_iter = 1,
                                    int max_iter = 200,
                                    double rho = 1,
                                    double alpha = 0.95,
                                    double trans_sd_init = 0.1) {
  int T = Y.size();
  int n = Y[0].n_rows;
  
  // Target parameters
  double mean_beta = 0;
  double sigma_beta = 0.5;
  double mean_beta_new = 0;
  double sigma_beta_new = 0.5;
  double trans_sd = trans_sd_init;
  double sigma_X1 = 0.2;
  
  // Initialize Mean_X as a vector of matrices (list of length T)
  std::vector<arma::mat> Mean_X(T);
  for (int t = 0; t < T; ++t) {
    Mean_X[t] = arma::mat(n, d, arma::fill::randn);
  }
  
  // Initialize Sigma_X as a vector of vectors of matrices
  std::vector<std::vector<arma::mat>> Sigma_X(T, std::vector<arma::mat>(n));
  for (int t = 0; t < T; ++t) {
    for (int i = 0; i < n; ++i) {
      Sigma_X[t][i] = arma::mat(d, d, arma::fill::eye);
    }
  }
  
  // Initialize Xi as a vector of matrices
  std::vector<arma::mat> Xi(T);
  for (int t = 0; t < T; ++t) {
    Xi[t] = arma::mat(n, n, arma::fill::ones);
  }
  
  int K = 500;
  arma::vec err(K, arma::fill::zeros);
  arma::vec train_auc(K, arma::fill::zeros);
  train_auc(0) = -10;
  
  int k = 1;
  int ind = 0;
  
  while (k < K && ind == 0) {
    double V_beta_cumulative = 0;
    double M_beta_cumulative = 0;
    
    // Initialize M_X_cumulative as a vector of matrices
    std::vector<arma::mat> M_X_cumulative(T);
    std::vector<std::vector<arma::mat>> V_X_cumulative(T, std::vector<arma::mat>(n));
    for (int t = 0; t < T; ++t) {
      M_X_cumulative[t] = arma::mat(n, d, arma::fill::zeros);
      for (int i = 0; i < n; ++i) {
        V_X_cumulative[t][i] = arma::mat(d, d, arma::fill::zeros);
      }
    }
    
    // Update Xi and auxiliary values
    for (int t = 0; t < T; ++t) {
      for (int i = 0; i < n; ++i) {
        arma::vec M_i = Mean_X[t].row(i).t();
        arma::mat V_i = Sigma_X[t][i];
        for (int j = 0; j < n; ++j) {
          arma::vec M_j = Mean_X[t].row(j).t();
          arma::mat V_j = Sigma_X[t][j];
          
          double c = std::sqrt(distance_squared_inner_prod(M_i, M_j, V_i, V_j) + mean_beta * mean_beta + sigma_beta * sigma_beta);
          Xi[t](i, j) = 1.0 / (2.0 * c) * (std::exp(c) - 1.0) / (std::exp(c) + 1.0) / (-2.0);
          
          if (j != i) {
            V_beta_cumulative -= 2.0 * Xi[t](i, j) * alpha;
            M_beta_cumulative += (Y[t](i, j) - 0.5 + 2.0 * Xi[t](i, j) * dot(M_i, M_j)) * alpha;
          }
        }
      }
    }
    
    // Update of sigma_beta
    sigma_beta_new = 1.0 / std::sqrt(1.0 / (sigma_beta_prior * sigma_beta_prior) + V_beta_cumulative);
    
    // Update of mean_beta
    mean_beta_new = sigma_beta_new * sigma_beta_new * (1.0 / (sigma_beta_prior * sigma_beta_prior) * mean_beta_prior + M_beta_cumulative);
    
    // Prepare updated Mean_X and Sigma_X
    std::vector<arma::mat> Mean_X_new(T);
    std::vector<std::vector<arma::mat>> Sigma_X_new(T, std::vector<arma::mat>(n));
    for (int t = 0; t < T; ++t) {
      Mean_X_new[t] = arma::mat(n, d, arma::fill::zeros);
      for (int i = 0; i < n; ++i) {
        Sigma_X_new[t][i] = arma::mat(d, d, arma::fill::zeros);
      }
    }
    
    // Initialize Forward and Backward variables
    std::vector<std::vector<arma::mat>> Forward_var(T, std::vector<arma::mat>(n));
    std::vector<std::vector<arma::mat>> Backward_var(T, std::vector<arma::mat>(n));
    std::vector<arma::mat> Forward_Mean(T);
    std::vector<arma::mat> Backward_Mean(T);
    for (int t = 0; t < T; ++t) {
      Forward_Mean[t] = arma::mat(n, d, arma::fill::zeros);
      Backward_Mean[t] = arma::mat(n, d, arma::fill::zeros);
      for (int i = 0; i < n; ++i) {
        Forward_var[t][i] = arma::mat(d, d, arma::fill::zeros);
        Backward_var[t][i] = arma::mat(d, d, arma::fill::zeros);
      }
    }
    
    std::vector<std::vector<arma::mat>> Cross_cov(T - 1, std::vector<arma::mat>(n, arma::mat(2 * d, 2 * d, arma::fill::zeros)));
    
    arma::mat I_d = arma::eye<arma::mat>(d, d); // Identity matrix of size d
    
    // Update for each node i
    for (int i = 0; i < n; ++i) {
      // Reset cumulative variables for node i
      for (int t = 0; t < T; ++t) {
        V_X_cumulative[t][i].zeros();
        M_X_cumulative[t].row(i).zeros();
      }
      
      // Update V_X_cumulative and M_X_cumulative
      for (int t = 0; t < T; ++t) {
        arma::vec M_i = Mean_X[t].row(i).t();
        arma::mat V_i = Sigma_X[t][i];
        for (int j = 0; j < n; ++j) {
          if (j != i) {
            arma::vec M_j;
            arma::mat V_j;
            if (j > i) {
              M_j = Mean_X[t].row(j).t();
              V_j = Sigma_X[t][j];
            } else {
              M_j = Mean_X_new[t].row(j).t();
              V_j = Sigma_X_new[t][j];
            }
            V_X_cumulative[t][i] -= 2.0 * Xi[t](i, j) * alpha * (M_j * M_j.t() + V_j);
            M_X_cumulative[t].row(i) += alpha * (Y[t](i, j) - 0.5 + 2.0 * Xi[t](i, j) * mean_beta_new) * M_j.t();
          }
        }
      }
      
      // Initialize the first forward variables
      Forward_var[0][i] = -std::pow(rho, -2) * std::pow(trans_sd, 4) *
        (std::pow(sigma_X1, -2) * arma::eye(d, d) + rho * rho * std::pow(trans_sd, -2) * arma::eye(d, d) + V_X_cumulative[0][i]);
      Forward_Mean[0].row(i) = (-rho * std::pow(trans_sd, 2)) * M_X_cumulative[0].row(i);
      
      // Initialize the last backward variables
      Backward_var[T - 1][i] = -std::pow(rho, -2) * std::pow(trans_sd, 4) *
        (std::pow(trans_sd, -2) * arma::eye(d, d) + V_X_cumulative[T - 1][i]);
      Backward_Mean[T - 1].row(i) = (-rho * std::pow(trans_sd, 2)) * M_X_cumulative[T - 1].row(i);
      
      // Forward and backward message passing
      if (T > 2) {
        for (int t2 = 1; t2 < T - 1; ++t2) {
          // Forward update
          Forward_var[t2][i] = -std::pow(rho, -2) * std::pow(trans_sd, 4) *
            (pinv(Forward_var[t2 - 1][i]) + (1.0 + rho * rho) * std::pow(trans_sd, -2) * arma::eye(d, d) + V_X_cumulative[t2][i]);
          Forward_Mean[t2].row(i) = (-rho * std::pow(trans_sd, 2)) *
            (pinv(Forward_var[t2 - 1][i]) * Forward_Mean[t2 - 1].row(i).t() + M_X_cumulative[t2].row(i).t()).t();
          
          // Backward update
          Backward_var[T - t2 - 1][i] = -std::pow(rho, -2) * std::pow(trans_sd, 4) *
            (pinv(Backward_var[T - t2][i]) + (1.0 + rho * rho) * std::pow(trans_sd, -2) * arma::eye(d, d) + V_X_cumulative[T - t2 - 1][i]);
          Backward_Mean[T - t2 - 1].row(i) = (-rho * std::pow(trans_sd, 2)) *
            (pinv(Backward_var[T - t2][i]) * Backward_Mean[T - t2].row(i).t() + M_X_cumulative[T - t2 - 1].row(i).t()).t();
        }
      }
      
      // Update Mean_X_new and Sigma_X_new
      for (int t = 0; t < T; ++t) {
        arma::mat temp_var;
        arma::vec temp_mean;
        if (t == 0) {
          temp_var = std::pow(sigma_X1, -2) * arma::eye(d, d) + rho * rho * std::pow(trans_sd, -2) * arma::eye(d, d) +
            pinv(Backward_var[t + 1][i]) + V_X_cumulative[t][i];
          temp_mean = pinv(Backward_var[t + 1][i]) * Backward_Mean[t + 1].row(i).t() + M_X_cumulative[t].row(i).t();
        } else if (t == T - 1) {
          temp_var = std::pow(trans_sd, -2) * arma::eye(d, d) + pinv(Forward_var[t - 1][i]) + V_X_cumulative[t][i];
          temp_mean = pinv(Forward_var[t - 1][i]) * Forward_Mean[t - 1].row(i).t() + M_X_cumulative[t].row(i).t();
        } else {
          temp_var = (1.0 + rho * rho) * std::pow(trans_sd, -2) * arma::eye(d, d) +
            pinv(Forward_var[t - 1][i]) + pinv(Backward_var[t + 1][i]) + V_X_cumulative[t][i];
          temp_mean = pinv(Forward_var[t - 1][i]) * Forward_Mean[t - 1].row(i).t() +
            pinv(Backward_var[t + 1][i]) * Backward_Mean[t + 1].row(i).t() + M_X_cumulative[t].row(i).t();
        }
        
        Sigma_X_new[t][i] = pinv(temp_var);
        Sigma_X_new[t][i] = 0.5 * (Sigma_X_new[t][i] + Sigma_X_new[t][i].t()); // Ensure symmetry
        Mean_X_new[t].row(i) = (Sigma_X_new[t][i] * temp_mean).t();
      }



      // First time point
      Cross_cov[0][i].submat(0, 0, d - 1, d - 1) = (1.0 / std::pow(sigma_X1, 2)) * I_d +
        (1.0 / std::pow(trans_sd, 2)) * I_d +
        V_X_cumulative[0][i];
      Cross_cov[0][i].submat(0, d, d - 1, 2 * d - 1) = -(1.0 / std::pow(trans_sd, 2)) * I_d;
      Cross_cov[0][i].submat(d, 0, 2 * d - 1, d - 1) = -(1.0 / std::pow(trans_sd, 2)) * I_d;
      Cross_cov[0][i].submat(d, d, 2 * d - 1, 2 * d - 1) = (1.0 / std::pow(trans_sd, 2)) * I_d +
        arma::inv(Backward_var[2][i]) +
        (1.0 / std::pow(trans_sd, 2)) * I_d +
        V_X_cumulative[1][i];
      Cross_cov[0][i] = arma::inv(Cross_cov[0][i]);
      
      // Intermediate time points
      for (int t = 1; t < T - 2; ++t) {
        Cross_cov[t][i].submat(0, 0, d - 1, d - 1) = (1.0 / std::pow(trans_sd, 2)) * I_d +
          arma::inv(Forward_var[t - 1][i]) +
          (1.0 / std::pow(trans_sd, 2)) * I_d +
          V_X_cumulative[t][i];
        Cross_cov[t][i].submat(0, d, d - 1, 2 * d - 1) = -(1.0 / std::pow(trans_sd, 2)) * I_d;
        Cross_cov[t][i].submat(d, 0, 2 * d - 1, d - 1) = -(1.0 / std::pow(trans_sd, 2)) * I_d;
        Cross_cov[t][i].submat(d, d, 2 * d - 1, 2 * d - 1) = (1.0 / std::pow(trans_sd, 2)) * I_d +
          arma::inv(Backward_var[t + 2][i]) +
          (1.0 / std::pow(trans_sd, 2)) * I_d +
          V_X_cumulative[t + 1][i];
        Cross_cov[t][i] = arma::inv(Cross_cov[t][i]);
      }
      
      // Last time point
      Cross_cov[T - 2][i].submat(0, 0, d - 1, d - 1) = (1.0 / std::pow(trans_sd, 2)) * I_d +
        arma::inv(Forward_var[T - 2][i]) +
        (1.0 / std::pow(trans_sd, 2)) * I_d +
        V_X_cumulative[T - 2][i];
      Cross_cov[T - 2][i].submat(0, d, d - 1, 2 * d - 1) = -(1.0 / std::pow(trans_sd, 2)) * I_d;
      Cross_cov[T - 2][i].submat(d, 0, 2 * d - 1, d - 1) = -(1.0 / std::pow(trans_sd, 2)) * I_d;
      Cross_cov[T - 2][i].submat(d, d, 2 * d - 1, 2 * d - 1) = (1.0 / std::pow(trans_sd, 2)) * I_d +
        V_X_cumulative[T - 1][i];
      Cross_cov[T - 2][i] = arma::inv(Cross_cov[T - 2][i]);
    
      }
    
    // Compute b and b0
    double b = 0.0;
    double b0 = 0.0;
    
    for (int i0 = 0; i0 < n; ++i0) {
      for (int t = 0; t < T - 1; ++t) {
        arma::rowvec Mean_t = (Mean_X_new[t]).row(i0);
        arma::rowvec Mean_t1 = (Mean_X_new[t + 1]).row(i0);
        arma::mat Cov_t = Cross_cov[t][i0];
        
        double EX = arma::accu(arma::square(Mean_t1)) + arma::accu(arma::square(Mean_t)) -
          2 * arma::dot(Mean_t, Mean_t1) +
          arma::trace(Cov_t.submat(0, 0, d - 1, d - 1) +
          Cov_t.submat(d, d, 2 * d - 1, 2 * d - 1) -
          2 * Cov_t.submat(0, d, d - 1, 2 * d - 1));
        
        b += EX;
      }
      
      arma::rowvec Mean_1 = (Mean_X_new[0]).row(i0);
      arma::mat Cov_1 = Sigma_X_new[0][i0];
      
      b0 += arma::accu(arma::square(Mean_1)) + arma::trace(Cov_1);
    }
    
    double tau_inv_sq = (n * (T - 1) * d + 1.0) / (b + 1.0);
    trans_sd = 1.0 / std::sqrt(tau_inv_sq);
    
    double sg_inv_sq = (n * d + 1.0) / (b0 + 1.0);
    sigma_X1 = 1.0 / std::sqrt(sg_inv_sq);
    
    // Compute error for convergence check
    arma::vec auc_mean((n - 1) * n * T / 2, arma::fill::zeros);
    arma::vec rmse_iter((n - 1) * n * T / 2, arma::fill::zeros);
    int r = 0;
    for (int t = 0; t < T; ++t) {
      for (int i = 0; i < n; ++i) {
        for (int j = 0; j < i; ++j) {
          double pred_new = 1.0 / (1.0 + std::exp(-mean_beta_new - arma::dot(Mean_X_new[t].row(i), Mean_X_new[t].row(j))));
          double pred_old = 1.0 / (1.0 + std::exp(-mean_beta - arma::dot(Mean_X[t].row(i), Mean_X[t].row(j))));
          auc_mean(r) = pred_new;
          rmse_iter(r) = pred_old;
          r++;
        }
      }
    }
    
    double temp_s = arma::dot(auc_mean - rmse_iter, auc_mean - rmse_iter);
    err(k - 1) = std::sqrt(temp_s / ((n - 1) * n * T / 2));
    
    if ((std::abs(err(k - 1)) < gap) || (k > max_iter)) {
      ind = 1;
    } else {
      printf("Iteration %d: Error = %f\n", k, err(k - 1));
      Mean_X = Mean_X_new;
      Sigma_X = Sigma_X_new;
      mean_beta = mean_beta_new;
      sigma_beta = sigma_beta_new;
      k++;
    }
  }
  
  // Prepare output
  Rcpp::List Mean_X_list(T);
  for (int t = 0; t < T; ++t) {
    Mean_X_list[t] = Mean_X[t];
  }
  Rcpp::List Sigma_X_list(T);
  for (int t = 0; t < T; ++t) {
    Rcpp::List Sigma_X_t_list(n);
    for (int i = 0; i < n; ++i) {
      Sigma_X_t_list[i] = Sigma_X[t][i];
    }
    Sigma_X_list[t] = Sigma_X_t_list;
  }
  
  return Rcpp::List::create(Named("Mean_X") = Mean_X_list,
                            Named("Sigma_X") = Sigma_X_list,
                            Named("iter") = k - 1,
                            Named("mean_beta") = mean_beta,
                            Named("sigma_beta") = sigma_beta,
                            Named("err") = err.head(k - 1));
}
