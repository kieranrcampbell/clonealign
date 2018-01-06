#include <Rcpp.h>
using namespace Rcpp;

using namespace std;

// [[Rcpp::export]]
double Qcpp(NumericVector pars,
         NumericMatrix gamma,
         NumericMatrix Y,
         NumericMatrix L,
         NumericVector S) {
  
  // Get sizes first
  int N = Y.nrow();
  int G = Y.ncol();
  int C = L.ncol();
  
  // Parse pars into a usable format
  NumericVector mu(G, 1.0);
  NumericVector phi(G, 1.0);
  
  for(int g = 0; g < G; g++) {
    if(g > 0)
      mu[g] = pars[g-1];
    phi[g] = pars[g + G - 1];
  }
  
  
  double qq = 0;
  
  NumericMatrix m2(G, C); // Mean of gene g on clone c
  
  for(int c = 0; c < C; c++) {
    NumericVector m(G, 0.0); // mean of clone c for g genes *unnormalised*
    for(int g = 0; g < G; g++) {
      m[g] = mu[g] * L(g,c);
    }
    
    // Compute sum
    double m_sum = 0;
    for(int g = 0; g < G; g++)
      m_sum += m[g];
    
    // Normalise
    for(int g = 0; g < G; g++)
      m2(g,c) = m[g] / m_sum;
  }
  
  // Now compute Q
  for(int c = 0; c < C; c++) {
    for(int n = 0; n < N; n++) {
      double l_nc = 0;
      for(int g = 0; g < G; g++)
        l_nc += gamma(n,c) * dnbinom_mu(Y(n,g), phi[g], m2(g,c) * S[n], 1);
      qq += l_nc;
    }
  }
  return -qq;
}


  