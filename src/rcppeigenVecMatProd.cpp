#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::RowVectorXd;
using Rcpp::as;
using Rcpp::wrap;
using Rcpp::List;
using Rcpp::NumericMatrix;
using Rcpp::NumericVector;

// [[Rcpp::export]]
List vectorMatrixProduct(NumericMatrix AA, NumericVector BB) {
	// Map the double matrix AA from R
	const Map<MatrixXd> A(as<Map<MatrixXd> >(AA));
	// Map the double vector BB from R
	const Map<RowVectorXd> B(as<Map<RowVectorXd> >(BB));
	return List::create(Rcpp::Named("result") = B*A);
}

