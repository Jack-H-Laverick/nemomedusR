#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// Index based weighted means from arrays in C++
// [[Rcpp::export]]
NumericVector array_w_mean(arma::cube array, DataFrame df) {
  
  int n = df.nrow();                           // Get df length to use as a counter
  
  NumericVector x = df["x"];                   // Make each column in the dataframe accessible
  NumericVector y = df["y"];
  NumericVector depth = df["depth"];
  NumericVector voxel = df["voxel"];
  NumericVector weight = df["weight"];
  
  int SummarySize = unique(voxel).size();      // Count the summaries to be produced
  
  NumericVector total(SummarySize);            // Start a vector to hold weighted totals per summary
  NumericVector total_w(SummarySize);          // Start a vector to hold total weights per summary
  
  for(int i = 0; i < n; ++i) {                 // From each row of the dataframe
    // I subtract 1 from each index because C++ starts at 0 (R starts at 1)
    total((voxel[i]-1)) += array((x[i]-1), (y[i]-1), (depth[i]-1)) * weight[i]; // Grab a value from the array and scale by weight, add this to the running total by summary
    total_w((voxel[i]-1)) += weight[i];                                         // Add the weight to the total weight for scaling the mean per summary
    
  }
  
  return total / total_w;                      // Return the weighted mean by dividing weighted total by total weights
  
}