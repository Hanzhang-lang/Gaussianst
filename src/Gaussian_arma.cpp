//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;
using namespace std;

//'@title test_func
//'@description
//'To check whether build package successfully
//'@param x  mat
//'@param n int
//'@import Rcpp
//'@import RcppArmadillo
//[[Rcpp::export()]]
bool vec1(arma::mat x, int n) {
  arma::mat A(1, 5, fill::randu);
  arma::mat B(5, 1, fill::randu);
  arma::mat C(5, 5, fill::randu);

  std::cout << A.is_vec() << endl;
  std::cout << B.is_vec() << endl;
  std::cout << C.is_vec() << endl;
  std::cout << x.is_vec() <<endl;
  return x.is_vec();
}

//'@title exponential_cov func
//'@description
//'calculate the exponential covariance function
//'@param x  vec
//'@param y vec
//'@useDynLib Gaussianst
//'@param params list for theta
//'@import Rcpp
//'@import RcppArmadillo
//'@export
//[[Rcpp::export()]]
arma::mat exponential_cov(arma::vec x, arma::vec y, List params){
  int n = y.size();
  std::cout << "n:"<< n <<endl;

  int m = x.size();
  arma::mat out(m, n);
  for(int i=0; i<m; i++){
    for (int j=0; j<n; j++){
      out(i, j) = double(params[0]) * exp( -0.5 * pow(x[i] - y[j], 2) * double(params[1]));
    }
  }
  return out;
    // double(params[0]) * exp( -0.5 * double(params[1]) * (rep(x, n).resize(m , n) - y)^2);
}

//'@title Gaussian_kernel func
//'@description
//'Calculate the Gaussian_kernel function using RcppArmadillo
//'@param x  mat X(m, p)
//'@param y mat Y(n, p)
//'@param kernel_width double
//'@useDynLib Gaussianst
//'@import Rcpp
//'@import RcppArmadillo
//'@examples
//'a = matrix(c(1, 2, 3, 4), nrow=2)
//'b = matrix(c(2, 3, 4, 5), nrow=2)
//'Gaussian_kernel(a, b, kernel_width = 1.0)
//'@export
//[[Rcpp::export()]]
arma::mat Gaussian_kernel(arma::mat x, arma::mat y, double kernel_width){
// ensure multivariate input: X(m,p), Y(n,p)

double s = kernel_width;
int m = size(x)[0];
int n = size(y)[0];
arma::mat out(m, n);
  for (int i=0; i<m; i++){
    for (int j=0; j<n; j++){
      out(i, j) = exp( -sum(pow(x.row(i) - y.row(j), 2)) / (s * s));
    }
  }
  return out;
}

