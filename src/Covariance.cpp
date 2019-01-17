// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

//' Correlation
//' 
//' Calculates the correlation between two vectors.
//' 
//' @param A First matrix.
//' @param B Second matrix.
//' @param cor Return correlation matrix?
//' @export
//' @return Numeric matrix. 
// [[Rcpp::export]]

SEXP cov(const Eigen::Map<Eigen::MatrixXd> A,const Eigen::Map<Eigen::MatrixXd> B,const bool cor = false){
  // Obs
  const int n = A.rows();
  // Cols
  const int p1 = A.cols();
  const int p2 = B.cols();
  // One vector
  const Eigen::VectorXd u = Eigen::VectorXd::Constant(n,1);
  // Center A
  Eigen::MatrixXd Za = Eigen::MatrixXd::Constant(n,p1,0);
  for(int j=0; j<p1; j++){
    Za.col(j) = A.col(j)-u*(A.col(j).mean());
  };
  // Center B
  Eigen::MatrixXd Zb = Eigen::MatrixXd::Constant(n,p2,0);
  for(int j=0; j<p2; j++){
    Zb.col(j) = B.col(j)-u*(B.col(j).mean());
  };
  // Scale
  if(cor){
    Za.colwise().normalize();
    Zb.colwise().normalize();
  };
  // Inner product
  Eigen::MatrixXd R = Eigen::MatrixXd::Constant(p1,p2,0);
  if(cor){
    R = Za.transpose()*Zb;
  } else {
    R = (Za.transpose()*Zb)/(n-1);
  };
  return Rcpp::wrap(R);
} 