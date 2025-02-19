# ifndef __ESTIMATE_PMR_BETA_STARTS_H
# define __ESTIMATE_PMR_BETA_STARTS_H


#include <RcppArmadillo.h>
#include <vector>
#include <math.h>
#include <algorithm>


using namespace Rcpp;


// Internal C++ code for use in functions exported to R
arma::vec beta_starts_cpp(const double& shape,
                          const double& offset,
                          const double& total0,
                          const uint32_t& compartments) {

    arma::vec times = arma::linspace<arma::vec>(0, 1, compartments+1U);
    times += offset;

    arma::vec stage_abunds0(compartments, arma::fill::none);
    double sum_stage_abunds0 = 0;

    double corr_time, pbeta_val, pbeta_val0;
    for (uint32_t i = 0; i <= compartments; i++) {
        corr_time = times(i);
        if (corr_time > 1) corr_time -= 1;
        pbeta_val = R::pbeta(corr_time, shape, shape, true, false);
        if (times(i) > 1) pbeta_val += 1;
        if (i > 0) {
            stage_abunds0(i-1) = total0 * (pbeta_val - pbeta_val0);
            sum_stage_abunds0 += stage_abunds0(i-1);
        }
        pbeta_val0 = pbeta_val;
    }

    double sum_diff = std::round(sum_stage_abunds0) != std::round(total0);
    if (sum_diff != 0) {
        stage_abunds0 *= (total0 / sum_stage_abunds0);
        sum_diff = std::round(arma::accu(stage_abunds0)) != total0;
    }
    if (sum_diff != 0) {
        std::string err = "beta_starts magnitude error (";
        err += std::to_string(sum_diff) + ", shape = ";
        err += std::to_string(shape) + ", offset = ";
        err += std::to_string(offset) + ")";
        Rcpp::warning(err.c_str());
    }
    return stage_abunds0;
}



#endif
