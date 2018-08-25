// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

//' Weighted Covariance
//'
//' Calculates the weighted covariance. 
//'
//' @param Y Numeric matrix.
//' @param m Mean vector.
//' @param w Weights.
//' @return Numeric matrix. 
// [[Rcpp::export]]
SEXP wCov(const Eigen::Map<Eigen::MatrixXd> Y, const Eigen::Map<Eigen::RowVectorXd> m, 
          const Eigen::Map<Eigen::VectorXd> w){
  // Dimensions
  const int n = Y.rows();
  const int p = Y.cols();
  // Denominator
  const double d = w.sum();
  // Numerator
  Eigen::MatrixXd A(p,p);
  A.setZero();
  Eigen::RowVectorXd z(p);
  z.setZero();
  
  for(int i=0; i<n; i++){
    z = (Y.row(i)-m);
    A += w(i)*(z.transpose()*z);
    // std::cout << A << std::endl;
  };
  // Result
  const Eigen::MatrixXd Out = (A/d);
  return Rcpp::wrap(Out);
}