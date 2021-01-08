#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

//' Index based weighted means from arrays in C++
//'
//' This function takes an array and a summary scheme and returns a vector of weighted means per group.
//' 
//' A summary scheme must contain 3 columns of array indices on rows, columns and slices, titled "x", "y", 
//' and "layer". A column named "group" defines the values to be summarised together. A column "weight" 
//' contains the weights for computing weighted mean. Helper functions are coming to help make common schemes.   
//'
//' @param array to be summarised.
//' @param scheme a dataframe containing the indices to average together with weights.
//' @return A vector containing summary values per group, ordered by group.
//' @export
// [[Rcpp::export]]
NumericVector array_w_mean(arma::cube array, DataFrame scheme) {
  
  int n = scheme.nrow();                       // Get scheme length to use as a counter
  
  NumericVector x = scheme["x"];               // Make each column in the dataframe accessible
  NumericVector y = scheme["y"];
  NumericVector layer = scheme["layer"];
  NumericVector group = scheme["group"];
  NumericVector weight = scheme["weight"];
  
  int SummarySize = unique(group).size();      // Count the summaries to be produced
  
  NumericVector total(SummarySize);            // Start a vector to hold weighted totals per summary
  NumericVector total_w(SummarySize);          // Start a vector to hold total weights per summary
  
  for(int i = 0; i < n; ++i) {                 // From each row of the dataframe
    // I subtract 1 from each index because C++ starts at 0 (R starts at 1)
    total((group[i]-1)) += array((x[i]-1), (y[i]-1), (layer[i]-1)) * weight[i]; // Grab a value from the array and scale by weight, add this to the running total by summary
    total_w((group[i]-1)) += weight[i];                                         // Add the weight to the total weight for scaling the mean per summary
    
  }
  
  return total / total_w;                      // Return the weighted mean by dividing weighted total by total weights
  
}