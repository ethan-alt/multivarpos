#include <RcppDist.h>
#include <RcppArmadillo.h>
#include <mvtnormAPI.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppDist)]]
// [[Rcpp::depends(mvtnorm)]]



// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
arma::vec kron_y(arma::mat const& Sigma, arma::mat const& Ymat) {
  int n  = Ymat.n_rows;
  int J  = Ymat.n_cols;
  arma::mat res(n, J, arma::fill::zeros);     // n x J matrix to store values
  for ( int j = 0; j < J; j++ ) {
    // arma::vec temp_j(n, arma::fill::zeros);   // temporary vector that will be jth column of res
    for ( int k = 0; k < J; k++ ) {
      // temp_j += Sigma[j,k] * Y.col(k);        // jth column = sum_{k=1}^J sigma[j,k] * Y[, k]
      res.col(j) += Sigma(j,k) * Ymat.col(k);
    }
    // res.col(j) = temp_j;
  }
  return arma::vectorise(res);
}

// [[Rcpp::export]]
arma::vec Xt_kron_y(arma::mat const& Sigma, List const& XtYlist, arma::vec const& pj) {
  int J = Sigma.n_cols;
  int p = sum(pj);
  arma::vec res(p, arma::fill::zeros);
  int start = 0;
  for ( int j = 0; j < J; j++ ) {
    int end = start + pj(j) - 1;
    arma::mat XtY_j = XtYlist[j];
    arma::vec temp_j(pj(j), arma::fill::zeros);
    for ( int k = 0; k < J; k++ ) {
      temp_j += Sigma(j,k) * XtY_j.col(k);
    }
    res.subvec(start, end) = temp_j;
    start = end + 1;
  }
  return res;
}

// [[Rcpp::export]]
arma::mat Xt_kron_X(arma::mat const& Sigma, arma::mat const& XtX, arma::vec const& pj) {
  int J = Sigma.n_cols;
  arma::mat res = XtX;
  int col_start = 0;
  for ( int j = 0; j < J; j++ ) {
    int row_start = 0;
    for ( int k = 0; k < J; k++ ) {
      res( row_start, col_start, size(pj(k), pj(j)) ) *= Sigma(j,k);
      row_start += pj(k);
    }
    col_start += pj(j);
  }
  return res;
}


// [[Rcpp::export]]
arma::mat get_A(arma::mat const& Y, arma::mat const& Xwide, arma::vec const& beta, arma::vec const& pj) {
  int n = Y.n_rows;
  int J = Y.n_cols;
  arma::mat A(n,J,arma::fill::zeros);   // initialize A = (Y - X*B), B = diag(beta_1, ..., beta_J)
  int first_indx = 0;
  for ( int j = 0; j < J; j++ ) {
    // subset beta
    int last_indx = first_indx + pj(j) - 1;
    A.col(j) = Y.col(j) - Xwide.cols(first_indx, last_indx) * beta.subvec(first_indx, last_indx);
    first_indx = last_indx + 1;
  }
  return (A.t() * A);  // cross product of (Y - X*B)
}


// [[Rcpp::export]]
arma::mat rwish_frac(double const& df, arma::mat const& S) {
  int J = S.n_rows;
  arma::mat A = arma::chol(S, "lower");
  arma::mat B(S.n_rows, S.n_cols, arma::fill::zeros);
  for ( int j = 0; j < J; j++ ) {
    B(j,j) = sqrt( R::rchisq(df - j) );
    for ( int k = j+1; k < J; k++ ) {
      B(k,j) = R::rnorm(0.0, 1.0);
    }
  }
  A = A * B;
  return A * A.t();
}


// [[Rcpp::export]]
List sur_sample_gibbs_rcpp(
    mat const& Xwide,
    mat const& Y,
    mat const& XtX,
    List const& XtYlist,
    int const& d0,
    vec beta,
    int const& nsmpl,
    vec const& pj,
    int const& burn = 0,
    int const& thin = 1,
    bool keep_sigma = true
) {
  int n  = Y.n_rows;
  int J  = Y.n_cols;
  int p  = sum(pj);
  
  // compute posterior degrees of freedom
  int post_df = n + d0 - J - 1;

  // create containers for beta and Sigma (if applicable)
  mat beta_sample(p, nsmpl, arma::fill::zeros);
  int nsmpl_sigma = 0;
  if ( keep_sigma ) 
    nsmpl_sigma = nsmpl;
  cube sigma_sample(J, J, nsmpl_sigma);
  mat Sigma(J,J,arma::fill::zeros);

  // start burn-in period
  for ( int i = 0; i < burn; i++ ) {
    // sample sigma from inverse-Wishart
      mat A = get_A(Y, Xwide, beta, pj);      // compute A = (Y - X * B)'(Y - X * B)
      Sigma = riwish(post_df, A);

    // sample beta from normal
       // compute X' Sigmainv \otimes I X
       mat Sigmainv  = inv_sympd(Sigma);
       mat cov_beta  = inv_sympd( Xt_kron_X(Sigmainv, XtX, pj) );
       vec mean_beta = cov_beta * Xt_kron_y(Sigmainv, XtYlist, pj);
       beta          = mvnrnd(mean_beta, cov_beta);
  }

  // burn-in period ended; start (thinned) sampling
  for ( int i = 0; i < nsmpl; i++ ) {
    for ( int k = 1; k <= thin; k++ ) {
      // sample sigma from inverse-Wishart
      mat A = get_A(Y, Xwide, beta, pj);      // compute A = (Y - X * B)'(Y - X * B)
      Sigma = riwish(post_df, A);
      
      // sample beta from normal
      // compute X' Sigmainv \otimes I X
      mat Sigmainv  = inv_sympd(Sigma);
      mat cov_beta  = inv_sympd( Xt_kron_X(Sigmainv, XtX, pj) );
      vec mean_beta = cov_beta * Xt_kron_y(Sigmainv, XtYlist, pj);
      beta          = mvnrnd(mean_beta, cov_beta);
    }
    // thin period ended; store sigma if applicable
    if ( keep_sigma )
      sigma_sample.slice(i) = Sigma;
    // store beta
    beta_sample.col(i) = beta;
  }
  beta_sample = beta_sample.t();
  if (keep_sigma) {
    return List::create(
      _["beta_sample"]  = beta_sample,
      _["sigma_sample"] = sigma_sample
    );
  }
  return List::create(
    _["beta_sample"] = beta_sample
  );
}






// [[Rcpp::export]]
List sur_sample_gibbs_pp_rcpp(
    mat const& Xwide,
    mat const& X0wide,
    mat const& Y,
    mat const& Y0,
    mat const& XtX,
    mat const& X0tX0,
    List const& XtYlist,
    List const& X0tY0list,
    int const& d0,
    vec beta,
    int const& nsmpl,
    vec const& pj,
    int const& burn = 0,
    int const& thin = 1,
    bool keep_sigma = true,
    double const& a0 = 1
) {
  int n  = Y.n_rows;
  int n0 = Y0.n_rows;
  int J  = Y.n_cols;
  int p  = sum(pj);
  
  // compute posterior degrees of freedom
  int post_df = n + a0 * n0 + d0 - J - 1;
  
  // create containers for beta and Sigma (if applicable)
  mat beta_sample(p, nsmpl, arma::fill::zeros);
  int nsmpl_sigma = 0;
  if ( keep_sigma ) 
    nsmpl_sigma = nsmpl;
  cube sigma_sample(J, J, nsmpl_sigma);
  mat Sigmainv(J, J, arma::fill::zeros);
  
  // start burn-in period
  for ( int i = 0; i < burn; i++ ) {
    // sample sigmainv Wishart
    mat A    = get_A(Y, Xwide, beta, pj) + a0 * get_A(Y0, X0wide, beta, pj);      // compute A = (Y - X * B)'(Y - X * B) + a0 * (Y0 - X0 * B)'(Y0 - X0 * B) 
    Sigmainv = rwish(post_df, arma::inv_sympd(A) );
    
    // sample beta from normal
    mat cov_beta  = inv_sympd( Xt_kron_X(Sigmainv, XtX, pj) + a0 * Xt_kron_X(Sigmainv, X0tX0, pj) );
    vec mean_beta = cov_beta * ( Xt_kron_y(Sigmainv, XtYlist, pj) + a0 * Xt_kron_y(Sigmainv, X0tY0list, pj) );
    beta          = mvnrnd(mean_beta, cov_beta);
  }
  
  // burn-in period ended; start (thinned) sampling
  for ( int i = 0; i < nsmpl; i++ ) {
    for ( int k = 1; k <= thin; k++ ) {
      // sample sigmainv from Wishart
      mat A    = get_A(Y, Xwide, beta, pj) + a0 * get_A(Y0, X0wide, beta, pj);      // compute A = (Y - X * B)'(Y - X * B) + a0 * (Y0 - X0 * B)'(Y0 - X0 * B) 
      Sigmainv = rwish(post_df, arma::inv_sympd(A) );
      
      // sample beta from normal
      mat cov_beta  = inv_sympd( Xt_kron_X(Sigmainv, XtX, pj) + a0 * Xt_kron_X(Sigmainv, X0tX0, pj) );
      vec mean_beta = cov_beta * ( Xt_kron_y(Sigmainv, XtYlist, pj) + a0 * Xt_kron_y(Sigmainv, X0tY0list, pj) );
      beta          = mvnrnd(mean_beta, cov_beta);
    }
    // thin period ended; store sigma if applicable
    if ( keep_sigma )
      sigma_sample.slice(i) = arma::inv_sympd(Sigmainv);
    // store beta
    beta_sample.col(i) = beta;
  }
  beta_sample = beta_sample.t();
  if (keep_sigma) {
    return List::create(
      _["beta_sample"]  = beta_sample,
      _["sigma_sample"] = sigma_sample
    );
  }
  return List::create(
    _["beta_sample"] = beta_sample
  );
}

// [[Rcpp::export]]
arma::mat rmatnorm(
    arma::mat const& M,
    arma::mat const& U_chol,
    arma::mat const& V_chol
) {
  arma::mat X(M.n_rows, M.n_cols, arma::fill::randn);
  return M + U_chol * X * V_chol;
}



// [[Rcpp::export]]
List mvlm_sample_rcpp (
    arma::mat const& Y,
    arma::mat const& X,
    int const& d0,
    int const& nsmpl,
    bool keep_sigma = true
) {
  int n = Y.n_rows;
  int J = Y.n_cols;
  int p = X.n_cols;
  int sigma_df = n + d0 - J - p - 1;
  arma::mat XtX_inv = inv_sympd(X.t() * X);
  arma::mat Bhat = XtX_inv * X.t() * Y;
  arma::mat A = (Y - X * Bhat);
  A = A.t() * A;
  
  arma::mat chol_XtX_inv = arma::chol(XtX_inv, "lower");
  
  // containers for samples
  int nsmpl_sigma = 0;
  if ( keep_sigma )
    nsmpl_sigma = nsmpl;
  arma::cube sigma_smpl(J, J, nsmpl_sigma);
  arma::cube beta_smpl(p, J, nsmpl);
  
  for ( int i = 0; i < nsmpl; i++ ) {
    // draw sigma
    arma::mat sigma = riwish(sigma_df, A);
    if (keep_sigma)
      sigma_smpl.slice(i) = sigma;
    
    // draw beta
    arma::mat beta = rmatnorm(Bhat, chol_XtX_inv, chol(sigma, "upper"));
    beta_smpl.slice(i) = beta;
  }
  if ( keep_sigma ) {
    return List::create(
      _["beta_sample"] = beta_smpl,
      _["sigma_sample"] = sigma_smpl
    );
  }
  else {
    return List::create(
      _["beta_sample"] = beta_smpl
    );
  }
}





// [[Rcpp::export]]
List mvlm_sample_pp_rcpp (
    arma::mat const& Y,
    arma::mat const& X,
    arma::mat const& Y0,
    arma::mat const& X0,
    int const& d0,
    int const& nsmpl,
    double const& a0,
    bool keep_sigma = true
) {
  int n0 = Y0.n_rows;
  int n  = Y.n_rows;
  int J  = Y.n_cols;
  int p  = X.n_cols;
  
  // helper quantities
  arma::mat XtX   = X.t() * X;
  arma::mat XtX0  = X0.t() * X0;
  arma::mat Bhat  = inv_sympd(XtX) * (X.t() * Y);
  arma::mat Bhat0 = inv_sympd(XtX0) * (X0.t() * Y0);
  arma::mat A     = Y - X * Bhat;
  arma::mat A0    = Y0 - X0 * Bhat0;
  A  = A.t() * A;
  A0 = A0.t() * A0;
  // posterior parameters
  arma::mat cov_B      = arma::inv_sympd(XtX + a0 * XtX0);
  arma::mat Lambda     = a0 * cov_B * XtX0;
  arma::mat mean_B     = Lambda * Bhat0 + (arma::eye(p,p) - Lambda) * Bhat;
  arma::mat chol_cov_B = arma::chol(cov_B, "lower");
  
  int sigma_df = n + a0 * n0 + d0 - J - p - 1;
  arma::mat Vn = A + a0 * A0 + (Bhat - Bhat0).t() * ( Lambda.t() * XtX ) * (Bhat - Bhat0);
  Rcout << "mean = \n" << Vn / (1.0 * (sigma_df - J - 1.0)) << "\n";
    
  // containers for samples
  int nsmpl_sigma = 0;
  if ( keep_sigma )
    nsmpl_sigma = nsmpl;
  arma::cube sigma_smpl(J, J, nsmpl_sigma);
  arma::cube beta_smpl(p, J, nsmpl);
  
  for ( int i = 0; i < nsmpl; i++ ) {
    // draw sigma
    arma::mat sigma = riwish(sigma_df, Vn);
    if ( keep_sigma )
      sigma_smpl.slice(i) = sigma;
    
    // draw beta
    arma::mat beta = rmatnorm(mean_B, chol_cov_B, chol(sigma, "upper"));
    beta_smpl.slice(i) = beta;
  }
  
  if ( keep_sigma ) {
    return List::create(
      _["beta_sample"] = beta_smpl,
      _["sigma_sample"] = sigma_smpl
    );
  }
  else {
    return List::create(
      _["beta_sample"] = beta_smpl
    );
  }
}



