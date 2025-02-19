
#include <RcppArmadillo.h>
#include <vector>
#include <math.h>
#include <algorithm>

#include "beta_starts.h"

using namespace Rcpp;



//' Produce a vector of starting abundances based on a symmetrical beta distribution
//'
//' @param shape Single numeric giving the Beta distribution shape parameter.
//' @param offset Single numeric giving the Beta distribution offset parameter.
//' @param total0 Single numeric giving the total abundance across all stages.
//' @param compartments Integer giving the number of compartments for
//'     stage structure.
//'
//' @export
//'
//[[Rcpp::export]]
NumericVector beta_starts(const double& shape,
                          const double& offset,
                          const double& total0,
                          const int& compartments) {

    if (shape <= 0) stop("shape <= 0");
    if (offset < 0) stop("offset < 0");
    if (offset > 1) stop("offset > 1");
    if (total0 <= 0) stop("total0 <= 0");
    if (compartments <= 0) stop("compartments <= 0");

    arma::vec stage_abunds0 = beta_starts_cpp(shape, offset, total0,
                                              static_cast<uint32_t>(compartments));
    NumericVector out(stage_abunds0.begin(), stage_abunds0.end());

    return out;
}


