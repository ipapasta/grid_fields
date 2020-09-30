#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace std;
using namespace Rcpp;
#include<iostream>
#include<cmath>
#include<time.h>

// [[Rcpp::export]]
arma::sp_mat hessian_diag_and_upper_triangle(double beta, int p, arma::mat L, arma::mat A, arma::mat expAx) {
  arma::sp_mat hessX(p, p);
  int i,j;
  arma::mat o(1,1);  
  for(i=0; i<p; i++){
    for(j=0; j<p; j++){
      o    = -exp(beta)*( ( ( (A.col(i)) % (A.col(j)) % L ).t()) * expAx );
      hessX(i,j) = o(0,0);
    }
  }  
  return hessX;
}


// experimentation
// -----------------------------
// arma::mat myfun(int n){
//   arma::mat B(n, n);
//   for(int i=0; i<n; i++)
//     {
//       for(int j=0; j<n; j++)
//         {
//           B(i,j)=i*j;
//         }
//     }
//   return B;
// }
