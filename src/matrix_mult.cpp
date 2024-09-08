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

// Batch multiple matrix multiplications for "pre" and "post" operations
// [[Rcpp::export]]
Rcpp::List batch_matrix_operations_rc(const arma::mat& asy_lin_rep_ps,
                                 const arma::mat& wols_eX, const arma::mat& XpX_inv,
                                 const arma::mat& wols_eX_treat, const arma::mat& XpX_inv_treat,
                                 const arma::colvec& M1, const arma::colvec& M2,
                                 const arma::colvec& M3, const arma::colvec& mom) {

  // Perform matrix-vector multiplication for inf.cont.ps
  arma::colvec inf_cont_ps = asy_lin_rep_ps * M2;

  // Perform matrix-matrix multiplication for asy.lin.rep.ols
  arma::mat asy_lin_rep_ols = wols_eX * XpX_inv;

  // Perform matrix-matrix multiplication for asy.lin.rep.ols.treat
  arma::mat asy_lin_rep_ols_treat = wols_eX_treat * XpX_inv_treat;

  // Perform matrix-vector multiplication for inf.treat.or
  arma::colvec inf_treat_or = asy_lin_rep_ols * M1;

  // Perform matrix-vector multiplication for inf.cont.or
  arma::colvec inf_cont_or = asy_lin_rep_ols * M3;

  // Perform matrix-vector multiplication for inf.or
  arma::colvec inf_or = (asy_lin_rep_ols_treat - asy_lin_rep_ols) * mom;

  // Return results in a list
  return Rcpp::List::create(Rcpp::Named("inf_cont_ps") = inf_cont_ps,
                            Rcpp::Named("inf_treat_or") = inf_treat_or,
                            Rcpp::Named("inf_cont_or") = inf_cont_or,
                            Rcpp::Named("inf_or") = inf_or);
}
