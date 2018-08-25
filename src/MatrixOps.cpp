// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

//' Matrix Trace
//'
//' Calculates the trace of a matrix \eqn{A}.
//'
//' @param A Numeric matrix.
//' @return Scalar. 
// [[Rcpp::export]]
SEXP tr(const Eigen::Map<Eigen::MatrixXd> A){
  const double t = A.diagonal().sum();
  return Rcpp::wrap(t);
}

//' Matrix Matrix Product
//'
//' Calculates the product \eqn{AB}. 
//'
//' @param A Numeric matrix.
//' @param B Numeric matrix.
//' @return Numeric matrix. 
// [[Rcpp::export]]
SEXP fastMMp(const Eigen::Map<Eigen::MatrixXd> A, const Eigen::Map<Eigen::MatrixXd> B){
  const Eigen::MatrixXd C = A*B;
  return Rcpp::wrap(C);
}

//' Matrix Transpose
//'
//' Constructs \eqn{A'} from \eqn{A}.
//'
//' @param A Numeric matrix.
//' @return Numeric matrix.
// [[Rcpp::export]]
SEXP fastT(const Eigen::Map<Eigen::MatrixXd> A){
  const Eigen::MatrixXd At = A.transpose();
  return Rcpp::wrap(At);
}

//' Matrix Inner Product
//'
//' Calculates the product \eqn{A'B}.
//'
//' @param A Numeric matrix.
//' @param B Numeric matrix.
//' @return Numeric matrix. 
// [[Rcpp::export]]
SEXP fastIP(const Eigen::Map<Eigen::MatrixXd> A, const Eigen::Map<Eigen::MatrixXd> B){
  const Eigen::MatrixXd AtB = (A.transpose() * B);
  return Rcpp::wrap(AtB);
}

//' Matrix Inverse
//' 
//' Calcualtes \eqn{A^{-1}}.
//'
//' @param A Numeric matrix.
//' @return Numeric matrix. 
// [[Rcpp::export]]
SEXP fastInv(const Eigen::Map<Eigen::MatrixXd> A){
  const Eigen::MatrixXd Ai = A.completeOrthogonalDecomposition().pseudoInverse();
  return Rcpp::wrap(Ai);
}

//' Matrix Determinant
//'
//' Calculates the determinant of \eqn{A}.
//'
//' @param A Numeric matrix.
//' @return Scalar. 
// [[Rcpp::export]]
SEXP fastDet(const Eigen::Map<Eigen::MatrixXd> A){
  const double d = A.determinant();
  return Rcpp::wrap(d);
}

//' Fast Outer Product
//' 
//' Calculates the outer product \eqn{xy'}.
//' 
//' @param x Numeric vector.
//' @param y Numeric vector.
//' @return Numeric matrix.
// [[Rcpp::export]]
SEXP fastOP(const Eigen::Map<Eigen::VectorXd> x, const Eigen::Map<Eigen::VectorXd> y){
  const Eigen::MatrixXd Q = x * y.transpose();
  return Rcpp::wrap(Q);
}

//' Quadratic Form
//' 
//' Calculates the quadratic form \eqn{X'AX}.
//' 
//' @param X Numeric matrix.
//' @param A Numeric matrix.
//' @return Numeric matrix.
// [[Rcpp::export]]
SEXP fastQF(const Eigen::Map<Eigen::MatrixXd> X, const Eigen::Map<Eigen::MatrixXd> A){
  const Eigen::MatrixXd Q = X.transpose()*A*X;
  return Rcpp::wrap(Q);
}

//' Schur complement
//'
//' Calculates the efficient information \eqn{I_{11}-I_{12}I_{22}^{-1}I_{21}}. 
//'
//' @param I11 Information of target parameter
//' @param I22 Information of nuisance parameter
//' @param I12 Cross information between target and nuisance parameters
//' @return Numeric matrix. 
//'
// [[Rcpp::export]]
SEXP SchurC(const Eigen::Map<Eigen::MatrixXd> I11, const Eigen::Map<Eigen::MatrixXd> I22,
            const Eigen::Map<Eigen::MatrixXd> I12){
  // Kernel matrix
  const Eigen::MatrixXd K = I11 - I12 * I22.ldlt().solve(I12.transpose());
  return Rcpp::wrap(K);
}