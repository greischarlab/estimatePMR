#include <RcppArmadillo.h>
#include <cmath>

#include "ode.h"

using namespace Rcpp;



struct ConstantPMRgammaN {

    double cycleLength; // Length of one full cycle
    double mu;          // Non-sequestered iRBC mortality rate
    double museq;       // Sequestered iRBC mortality rate
    double R;           // Parasite multiplication rate
    uint32_t n;         // Number of age classes
    double lambdaN;     // Transition rate for non-sequestered iRBCs
    double lambdaS;     // Transition rate for sequestered iRBCs
    arma::vec qValues;

    ConstantPMRgammaN(const VecType& parms)
        : cycleLength(parms[0]),
          mu(parms[1]),
          museq(parms[2]),
          R(parms[3]),
          n(static_cast<uint32_t>(parms[4])),
          lambdaN(parms[4] / cycleLength),
          lambdaS(lambdaN),
          qValues(n, arma::fill::zeros) {

        double dt = 1 / lambdaN;
        double yValue0 = yfx(dt);
        double yValue;

        for (uint32_t i = 1; i <= n; ++i) {
            /*
             Update qValue.
             Note that this will naturally result in qValue = 0 when i = 1, which
             is what we want bc (1) there is no sequestration for the first
             age class and (2) for both dN and dS, we're using the qValue for
             the previous age class.
             */
            yValue = yfx(dt + (i-1) * dt);
            qValues(i-1) = 1 - yValue / yValue0; // Proportional change in yValues.
            yValue0 = yValue;

        }

    }


    void operator()(const VecType& x, VecType& dxdt, const double t) {

        // For the first age class (i = 0):
        dxdt[0] = R * (lambdaN * x[n - 1] + lambdaS * x[2 * n - 1]) -
            (lambdaN + mu) * x[0];
        dxdt[n] = -(lambdaS + museq) * x[n];

        for (uint32_t i = 1; i < n; ++i) {

            const double& Ni(x[i]);     // non-sequestered, age class i
            const double& Ni1(x[i-1]);  // non-sequestered, age class i-1
            double& dNi(dxdt[i]);       // change in non-sequestered, age class i

            const double& Si(x[n+i]);   // sequestered, age class i
            const double& Si1(x[n+i-1]);// sequestered, age class i-1
            double& dSi(dxdt[n+i]);     // change in sequestered, age class i
            const double& q(qValues(i-1));

            dNi = (1 - q) * lambdaN * Ni1 - (lambdaN + mu) * Ni;
            dSi = q * lambdaN * Ni1 + lambdaS * Si1 - (lambdaS + museq) * Si;

        }

        return;

    }


private:

    // Define the values of all parameters used inside yfx.
    const double p1 = 11.3869 / 467.6209; // Lower bound
    const double p2 = 1.0;                // Upper bound
    const double p3 = 18.5802 / 2.0;      // Inflection point
    const double p4 = 0.2242 * 2.0;       // Slope

    // gfx function: Computes the level of sequester protein expression at a certain iRBC age.
    // NOTE: yfx(age) = 1 - gfx(age)
    // The function uses a sigmoid curve defined by the given parameters.
    // Parameters p1, p2, p3, and p4 determine the shape and position of the curve.
    double yfx(const double& age) const {

      double yVal = 1 - (p1 + ((p2 - p1) / (1 + std::pow(10, p4 * (p3 - age)))));

      return yVal;
    }

};





arma::mat constPMR_gammaN_ode_cpp(std::vector<double> x0,
                                  const std::vector<double>& parms,
                                  const double& max_t,
                                  const double& dt) {

    if (max_t <= 0) stop("max_t must be > 0");
    if (dt >= max_t || dt <= 0) stop("dt must be < max_t and > 0");

    if (parms.size() != 5U) stop("`parms` must be of length 5");
    // The check below is important bc turning a negative double to an unsigned
    // integer will result in a very large number (in this case, over 4e9).
    if (parms[4] < 0) stop("5th element of `parms` cannot be < 0");

    uint32_t n = static_cast<uint32_t>(parms[4]);
    if (x0.size() != (2*n))
        stop("`x0` must be twice length as the 5th element of `parms`");

    GreedyObserver<VecType> obs;
    ConstantPMRgammaN system(parms);

    boost::numeric::odeint::integrate_const(
        VecStepperType(), std::ref(system),
        x0, 0.0, max_t, dt, std::ref(obs));

    uint32_t n_steps = obs.data.size();
    arma::mat output(n_steps, x0.size()+1U);
    for (uint32_t i = 0; i < n_steps; i++) {
        output(i,0) = obs.time[i];
        for (uint32_t j = 0; j < x0.size(); j++) {
            output(i,j+1U) = obs.data[i][j];
        }
    }

    return output;

}




//' Run ODE for constant PMR and gamma N.
//'
//' @param x0 Initial conditions for all stages (length should be `2*parms[5]`).
//' @param parms Vector of length 5 containing parameter values.
//'     Parameters in order: `cycleLength`, `mu`, `museq`, `R`, and `n`.
//' @param Maximum time point to simulate to.
//' @param Time step to use for ODE.
//'
//'
//' @export
//'
//' @return Numeric matrix where rows are time and columns are stages.
//'
//[[Rcpp::export]]
arma::mat constPMR_gammaN_ode(const std::vector<double>& x0,
                              const std::vector<double>& parms,
                              const double& max_t,
                              const double& dt) {

    arma::mat out = constPMR_gammaN_ode_cpp(x0, parms, max_t, dt);

    return out;

}


