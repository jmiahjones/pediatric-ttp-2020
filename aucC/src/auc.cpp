#include <Rcpp.h>
using namespace Rcpp;

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
double aucC(NumericVector pos_predictor, NumericVector neg_predictor){
  int pos_n = pos_predictor.size();
  int neg_n = neg_predictor.size();
  
  double sgn_total = 0;
  for(int pos_idx = 0; pos_idx < pos_n; ++pos_idx){
    for(int neg_idx = 0; neg_idx < neg_n; ++neg_idx){
      sgn_total += (1.0 + ((pos_predictor[pos_idx] - neg_predictor[neg_idx]) > 0) - ((pos_predictor[pos_idx] - neg_predictor[neg_idx]) < 0)) / 2.0;
    }
  }
  
  double auc = sgn_total / (pos_n * neg_n);
  return auc;
}


