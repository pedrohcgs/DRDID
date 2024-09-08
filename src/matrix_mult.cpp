// Include Rcpp and RcppArmadillo headers
#include <RcppArmadillo.h>

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

// Optimize crossproduct using RcppArmadillo
// [[Rcpp::export]]
arma::mat opt_crossprod(const arma::mat& A, const arma::mat& B, double n) {
  return (A.t() * B) / n;
}


// [[Rcpp::depends(RcppArmadillo)]]

// Batch multiple matrix multiplications
// [[Rcpp::export]]
Rcpp::List batch_matrix_operations(const arma::mat& wols_eX, const arma::mat& XpX_inv,
                                   const arma::mat& score_ps, const arma::mat& Hessian_ps,
                                   const arma::colvec& M1, const arma::colvec& M2, const arma::colvec& M3) {

  // First, perform matrix multiplications
  arma::mat asy_lin_rep_wols_result = wols_eX * XpX_inv;
  arma::mat asy_lin_rep_ps_result = score_ps * Hessian_ps;

  // Perform matrix-vector multiplications
  arma::colvec inf_treat_2_result = asy_lin_rep_wols_result * M1;
  arma::colvec inf_cont_2_result = asy_lin_rep_ps_result * M2;
  arma::colvec inf_cont_3_result = asy_lin_rep_wols_result * M3;

  // Return results in a list
  return Rcpp::List::create(Rcpp::Named("inf_treat_2") = inf_treat_2_result,
                            Rcpp::Named("inf_cont_2") = inf_cont_2_result,
                            Rcpp::Named("inf_cont_3") = inf_cont_3_result);
}
