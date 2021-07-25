// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

//' Matrix Determinant
//'
//' Calculates the determinant of \eqn{A}.
//'
//' @param A Numeric matrix.
//' @param logDet Return the logarithm of the determinant? 
//' @return Scalar. 
//' @noRd
// [[Rcpp::export]]
SEXP matDet(const arma::mat A, const bool logDet=false){
  double d;
  double s;
  if(logDet){
    arma::log_det(d, s, A);
  } else {
    d = arma::det(A);
  }
  return Rcpp::wrap(d);
}

//' Matrix Inverse
//' 
//' Calcualtes \eqn{A^{-1}}.
//'
//' @param A Numeric matrix.
//' @return Numeric matrix. 
//' @noRd
// [[Rcpp::export]]
SEXP matInv(const arma::mat A){
  const arma::mat Ai = arma::pinv(A);
  return Rcpp::wrap(Ai);
}

//' Matrix Inner Product
//'
//' Calculates the product \eqn{A'B}.
//'
//' @param A Numeric matrix.
//' @param B Numeric matrix.
//' @return Numeric matrix.
//' @noRd
// [[Rcpp::export]]
SEXP matIP(const arma::mat A, const arma::mat B){
  const arma::mat AtB = A.t() * B;
  return Rcpp::wrap(AtB);
}

//' Matrix Matrix Product
//'
//' Calculates the product \eqn{AB}. 
//'
//' @param A Numeric matrix.
//' @param B Numeric matrix.
//' @return Numeric matrix.
//' @noRd
// [[Rcpp::export]]
SEXP MMP(const arma::mat A, const arma::mat B){
  const arma::mat C = A * B;
  return Rcpp::wrap(C);
}

//' Matrix Outer Product
//' 
//' Calculates the outer product \eqn{AB'}.
//' 
//' @param A Numeric matrix.
//' @param B Numeric matrix.
//' @return Numeric matrix.
//' @noRd
// [[Rcpp::export]]
SEXP matOP(const arma::mat A, const arma::mat B){
  const arma::mat ABt = A * B.t();
  return Rcpp::wrap(ABt);
}

//' Quadratic Form
//' 
//' Calculates the quadratic form \eqn{X'AX}.
//' 
//' @param X Numeric matrix.
//' @param A Numeric matrix.
//' @return Numeric matrix.
//' @noRd
// [[Rcpp::export]]
SEXP matQF(const arma::mat X, const arma::mat A){
  const arma::mat xAx = X.t() * A * X;
  return Rcpp::wrap(xAx);
}

//' Schur complement
//'
//' Calculates the efficient information \eqn{I_{bb}-I_{ba}I_{aa}^{-1}I_{ab}}. 
//'
//' @param Ibb Information of target parameter
//' @param Iaa Information of nuisance parameter
//' @param Iba Cross information between target and nuisance parameters
//' @return Numeric matrix. 
//' @noRd
// [[Rcpp::export]]
SEXP SchurC(const arma::mat Ibb, const arma::mat Iaa,
            const arma::mat Iba){
  const arma::mat Ibba = Ibb - Iba * arma::solve(Iaa, Iba.t(), arma::solve_opts::likely_sympd);
  return Rcpp::wrap(Ibba);
}

//' Matrix Trace
//'
//' Calculates the trace of a matrix \eqn{A}.
//'
//' @param A Numeric matrix.
//' @return Scalar.
//' @noRd
// [[Rcpp::export]]
SEXP tr(const arma::mat A){
  const double t = arma::trace(A);
  return Rcpp::wrap(t);
}