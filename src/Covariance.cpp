// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

//' Covariance
//' 
//' Calculates the correlation between two matrices.
//' 
//' @param A NxP matrix.
//' @param B NxQ matrix.
//' @param corMat Return correlation matrix? If false, returns a covariance matrix.
//' @param eps Optional ridge parameter added to the diagonal of the covariance matrix.
//' @return Numeric matrix. 
//' @noRd
// [[Rcpp::export]]

SEXP matCov(
    const arma::mat A, 
    const arma::mat B, 
    const bool corMat=false,
    const double eps=0.0
){
  
  // Dimensions.
  const int n = A.n_rows;
  const int p = A.n_cols;
  const int q = B.n_cols;
  
  // Initialize.
  const arma::colvec u = arma::ones(n);
  arma::mat Za = arma::zeros(n, p);
  arma::mat Ma = arma::zeros(n, p);
  arma::mat Zb = arma::zeros(n, q);
  arma::mat Mb = arma::zeros(n, q);
  arma::mat R = arma::zeros(p, q);
  
  // Column-wise means.
  arma::rowvec ma = arma::mean(A, 0);
  arma::rowvec mb = arma::mean(B, 0);
  
  // Generate centering matrices.
  for(int i=0; i<n; i++){
  	Ma.row(i) = ma;
  	Mb.row(i) = mb;
  }
  
  // Centered matrices.
  Za = A - Ma;
  Zb = B - Mb;
  
  // Scale, if correlation matrix. 
  if(corMat){
    Za = arma::normalise(Za);
    Zb = arma::normalise(Zb);
  } 
  
  // Calculate covariance.
  R = Za.t() * Zb;
  if (not corMat) {
    R += eps * arma::eye(p, q);
    R /= (n-1);
  }
  
  // Output
  return Rcpp::wrap(R);
} 


//' Eigenvalues of Symmetric Matrix. 
//' 
//' Calculates the eigenvalues of a symmetric matrix. 
//' 
//' @param A symmetric matrix. 
//' @return Numeric vector.
//' @noRd
// [[Rcpp::export]]

SEXP eigSym(const arma::mat A) {
  
  arma::vec e;
  e = arma::eig_sym(A);
  
  // Output.
  return Rcpp::wrap(e);
} 