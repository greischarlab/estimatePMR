#include <RcppArmadillo.h>
// #include <cmath>
// #include <vector>
// #include <algorithm>
// #include <stdexcept>
// #include <iostream>
// #include <boost/math/distributions/beta.hpp>
// #include "ode.h"
// #include "beta_starts.h"
#include "archer_shared.h"


using namespace Rcpp;
using namespace arma;






//' Main optimization function for 8-parameter model
//'
//' @inheritParams archer_fit8_odeint
//'
//' @export
//'
// [[Rcpp::export]]
double archer_fit7_odeint(NumericVector parms, DataFrame data) {

  int n = get_n(parms[4]);
  //  Rcpp::Rcout << "n: " << n << std::endl;

  double pfCycleLength = 24;
  // double ageClassDuration = pfCycleLength / n;

  double varBetaDist = std::exp(-std::exp(parms[0])) / 12.0;

  // Rcpp::Rcout << "parms[0]: " << parms[0] << std::endl;
  // Rcpp::Rcout << "varBetaDist1: " << varBetaDist1 << std::endl;
  // Rcpp::Rcout << "varBetaDist: " << varBetaDist << std::endl;
  // Ensure varBetaDist is within a valid range
  //if (varBetaDist <= 0 || varBetaDist > 1) {
  //  stop("Invalid variance for Beta distribution: varBetaDist must be in (0, 1]");
  //}

  double betaShape = ((1 / varBetaDist) - 4) / 8;
  // Rcpp::Rcout << "betaShape: " << betaShape << std::endl;

  // Rcpp::Rcout << "betaShape: " << betaShape << std::endl;
  // Ensure betaShape is positive and finite
  // if (betaShape <= 0 || !std::isfinite(betaShape)) {
  //   stop("Invalid betaShape: must be positive and finite");
  // }


  double offset = std::exp(-std::exp(parms[1]));
  //  Rcpp::Rcout << "offset: " << offset << std::endl;

  double R = std::exp(parms[2]);
  // Rcpp::Rcout << "R: " << R << std::endl;

  double I0 = std::exp(parms[3]);
  //Rcpp::Rcout << "I0: " << I0 << std::endl;

  double inflec = std::exp(-std::exp(parms[5])) * 8 + 14;

  double ring_duration = std::exp(-std::exp(parms[6])) * 6 + 3;


  double sse = last_stage(pfCycleLength, R, n, inflec, betaShape, offset,
                          I0, ring_duration, data);

  return sse;


}
