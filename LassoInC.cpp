#include <cmath> // for math functions on scalars such as std::max and std::abs
#include <algorithm> // for std::max
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// Soft-thresholding function, returns scalar
// [[Rcpp::export]]
double soft_c(double a, double lambda){
  double sgn_a = (a > 0) - (a < 0); // returns +1 if a > 0, -1 if a < 0, and 0 if a == 0
  // Apply soft-thresholding: sign(a) * max(|a| - λ, 0)
  return sgn_a * std::max(std::abs(a) - lambda, 0.0);
}

// Lasso objective function, returns scalar
// [[Rcpp::export]]
double lasso_c(const arma::mat& Xtilde, const arma::colvec& Ytilde, const arma::colvec& beta, double lambda){
  // Number of observations
  double n = Xtilde.n_rows;
  // Compute residuals: Y - Xβ
  arma::colvec residuals = Ytilde - Xtilde * beta;
  // Compute residual sum of squares: ||Y - Xβ||²
  double rss = arma::dot(residuals, residuals);
  // Compute L1 penalty: sum of absolute coefficients
  double l1 = arma::norm(beta, 1);
  // Return LASSO objective: (1 / (2n)) * RSS + λ * ||β||₁
  return(rss / (2.0 * n) + lambda * l1);
}

// Lasso coordinate-descent on standardized data with one lamdba. Returns a vector beta.
// [[Rcpp::export]] 
arma::colvec fitLASSOstandardized_c(const arma::mat& Xtilde, const arma::colvec& Ytilde, double lambda, const arma::colvec& beta_start, double eps = 0.001){
  // Calculate initial objective value
  double f = lasso_c(Xtilde, Ytilde, beta_start, lambda); 
  // Initialize value of beta
  arma::colvec beta = beta_start;
  // Calculate number of rows and columns in Xtilde
  double n = Xtilde.n_rows;
  double p = Xtilde.n_cols;
  while (true) {
    // Update previous objective value
    double f_old = f;
    // Update previous value of beta
    arma::colvec beta_old = beta;
    
    arma::colvec XB = Xtilde * beta_old;
    // Calculate residual vector
    arma::colvec r = Ytilde - XB;
    for(int j = 0; j < p; j++) {
      // gradient component for coordinate j: (1/n) * X_jᵀR 
      double XjR_scaled = arma::dot(Xtilde.col(j), r) / n;
      // Update β_j via soft-thresholding: S(β_old[j] + (X_jᵀ r)/n, λ)
      beta(j) = soft_c(beta_old(j) + XjR_scaled, lambda);
      // Compute coefficient change for feature j: β_old[j] − β_new[j]
      double delta_beta_j = beta_old(j) - beta(j);
      // Update residuals: r ← r + X_j * (β_old[j] − β[j])
      r += Xtilde.col(j) * delta_beta_j;
    }
    // Evaluate LASSO objective at the updated coefficients β
    f = lasso_c(Xtilde, Ytilde, beta, lambda);
    // If the objective function doesn't decrease significantly, stop 
    if (f_old - f < eps) {
      break;
    }
  }
  return beta;
}  

// Lasso coordinate-descent on standardized data with supplied lambda_seq. 
// You can assume that the supplied lambda_seq is already sorted from largest to smallest, and has no negative values.
// Returns a matrix beta (p by number of lambdas in the sequence)
// [[Rcpp::export]]
arma::mat fitLASSOstandardized_seq_c(const arma::mat& Xtilde, const arma::colvec& Ytilde, const arma::colvec& lambda_seq, double eps = 0.001){
  // Retrieve dimensions of Xtilde
  int n = Xtilde.n_rows;
  int p = Xtilde.n_cols;
    
  // calculate length of lambda_seq
  int num_lambda = lambda_seq.n_elem;
  // initialize matrix of solutions at each lambda value
  arma::mat beta_mat(p, num_lambda);
  // Apply fitLASSOstandardized going from largest to smallest lambda 
  // (make sure supplied eps is carried over). 
  // Use warm starts strategy discussed in class for setting the starting values.
  arma::colvec beta = arma::colvec(p, arma::fill::zeros); 
  for (int i = 0; i < num_lambda; i++){
    double lambda = lambda_seq(i);
    // Perform coordinate descent LASSO for the current value of λ
    beta = fitLASSOstandardized_c(Xtilde, Ytilde, lambda, beta, eps);
    // Store the optimal coefficients and objective value for this λ
    beta_mat.col(i) =  beta;
  }
  return beta_mat;
}