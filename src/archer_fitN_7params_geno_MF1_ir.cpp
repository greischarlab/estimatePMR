#include <RcppArmadillo.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <boost/math/distributions/beta.hpp>
#include <ode.h>

using namespace Rcpp;
using namespace arma;

// Define the values of all parameters used inside yfx.
constexpr double p1 = 11.3869 / 467.6209; // Lower bound
constexpr double p2 = 1.0;               // Upper bound
constexpr double p4 = 0.2242 * 2.0;      // Slope

// Function to calculate yfx
NumericVector yfx(NumericVector age, double inflec) {
  int n = age.size();
  NumericVector yVals(n);
  
  for (int i = 0; i < n; ++i) {
    // Use Rcpp's `pow` for element-wise calculation
    yVals[i] = 1 - (p1 + ((p2 - p1) / (1 + pow(10, p4 * (inflec - age[i])))));
  }
  return yVals;
}

// Function to subset rows of a list
// Function to subset rows of a matrix and return a NumericMatrix
NumericMatrix subsetRows(const arma::mat& input, int step) {
  int n_rows = input.n_rows;
  if (n_rows == 0) {
    return NumericMatrix(0, 0);  // Return an empty matrix if input is empty
  } 
  
  // Calculate the number of rows in the subset
  int subset_n_rows = (n_rows -0 + step -1) / step;
  
  // Initialize a NumericMatrix to store the subset
  NumericMatrix subset(subset_n_rows, input.n_cols);
  
  // Fill the subset matrix
  for (int col = 0; col < input.n_cols; ++col) {
    arma::vec col_data = input.col(col);
    int row = 0; // Reset row counter for each column
    for (int i = 0; i < n_rows; i += step) {
      if (row >= subset_n_rows) {
        break; // Prevent out-of-bounds access
      }
      subset(row, col) = col_data(i);
      ++row;
    }
  }
  return subset;
} 


// beta_starts function betaShape, offset, I0, n
arma::vec beta_starts(const double& shape,
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


// constantPMR_gammaN function
struct ConstantPMRgammaN {
  
  double cycleLength; // Length of one full cycle
  double mu;          // Non-sequestered iRBC mortality rate
  double museq;       // Sequestered iRBC mortality rate
  double R;           // Parasite multiplication rate
  uint32_t n;         // Number of age classes
  double lambdaN;     // Transition rate for non-sequestered iRBCs
  double lambdaS;     // Transition rate for sequestered iRBCs
  double inflec;
  arma::vec qValues;
  
  ConstantPMRgammaN(const VecType& parms)
    : cycleLength(parms[0]),
      mu(parms[1]),
      museq(parms[2]),
      R(parms[3]),
      n(static_cast<uint32_t>(parms[4])),
      lambdaN(parms[4] / cycleLength),
      lambdaS(lambdaN),
      inflec(parms[5]),
      qValues(n, arma::fill::zeros) {
    
    double dt = 1 / lambdaN;
    double yValue0 = yfx1(dt, inflec);
    double yValue;
    
    for (uint32_t i = 1; i <= n; ++i) {
      yValue = yfx1(dt + (i-1) * dt, inflec);
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
  const double p1 = 11.3869/467.6209; // Lower bound
  const double p2 = 1.0;               // Upper bound
  const double p4 = 0.2242 * 2.0;       // Slope
  
  double yfx1(const double& age, double inflec) const {
    
    double yVal = 1 - (p1 + ((p2 - p1) / (1 + std::pow(10, p4 * (inflec - age)))));
    
    return yVal;
  }
  
};


arma::mat constPMR_gammaN_ode_cpp(std::vector<double> x0,
                                  const std::vector<double>& parms,
                                  const double& max_t,
                                  const double& dt) {
  
  if (max_t <= 0) stop("max_t must be > 0");
  if (dt >= max_t || dt <= 0) stop("dt must be < max_t and > 0");
  
  //if (parms.size() != 5U) stop("`parms` must be of length 5");
  // The check below is important bc turning a negative double to an unsigned
  // integer will result in a very large number (in this case, over 4e9).
  // if (parms[4] < 0) stop("5th element of `parms` cannot be < 0");
  
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


NumericVector repeat_subvector(const NumericVector& x) {
  std::vector<int> reps = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
                           6, 6, 6, 6, 6, 6, 6,
                           6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6};
  
  int totalSize = 146;
  
  // Create the output vector
  NumericVector result(totalSize);
  
  // Fill the output vector with repeated values
  int index = 0;
  for (int i = 0; i < x.size(); i++) {
    for (int j = 0; j < reps[i]; j++) {
      result[index++] = x[i];
    }
  } 
  
  return result;
}  


// Main optimization function
// [[Rcpp::export]] 
double archer_fitN_odeint(List parms, DataFrame data) {
  
  double cv_cycleLength = exp(-exp(as<double>(parms[4])));
  double squared = cv_cycleLength * cv_cycleLength;
  
  double reciprocal = 1.0 / squared;
  
  double result = std::round(reciprocal);
  
  int n = static_cast<int>(result);
  
  n = std::min(std::max(n, 4), 384);
  //  Rcpp::Rcout << "n: " << n << std::endl;
  
  double pfCycleLength = 24;
  double ageClassDuration = pfCycleLength / n;
  
  double varBetaDist1 = exp(-exp(as<double>(parms[0])));
  double varBetaDist = varBetaDist1 / 12;
  
  // Rcpp::Rcout << "parms[0]: " << parms[0] << std::endl;
  // Rcpp::Rcout << "varBetaDist1: " << varBetaDist1 << std::endl;
  // Rcpp::Rcout << "varBetaDist: " << varBetaDist << std::endl;
  // Ensure varBetaDist is within a valid range
  //if (varBetaDist <= 0 || varBetaDist > 1) {
  //  stop("Invalid variance for Beta distribution: varBetaDist must be in (0, 1]");
  //}
  
  double betaShape1 =  (1 / varBetaDist) - 4; 
  double betaShape = betaShape1/8;
  // Rcpp::Rcout << "betaShape: " << betaShape << std::endl;
  
  // Rcpp::Rcout << "betaShape: " << betaShape << std::endl;
  // Ensure betaShape is positive and finite
  // if (betaShape <= 0 || !std::isfinite(betaShape)) {
  //   stop("Invalid betaShape: must be positive and finite");
  // }
  
  
  double offset = exp(-exp(as<double>(parms[1])));
  //  Rcpp::Rcout << "offset: " << offset << std::endl;
  
  double R = exp(as<double>(parms[2]));
  // Rcpp::Rcout << "R: " << R << std::endl;
  
  double I0 = exp(as<double>(parms[3]));
  //Rcpp::Rcout << "I0: " << I0 << std::endl;
  
  double inflec = std::exp(-std::exp(as<double>(parms[5]))) * 8 + 14;
  
  double ring_duration = std::exp(-std::exp(as<double>(parms[6]))) * 6 + 3;
  
  // Initializing parameters
  std::vector<double> parmValues = {pfCycleLength, 0.0, 0.0, R, static_cast<double>(n), inflec};
  //  std::cout << "parmValues: ";
  //for (size_t i = 0; i < parmValues.size(); ++i) {
  //  std::cout << parmValues[i];
  //   if (i < parmValues.size() - 1) {
  //     std::cout << ", ";  // Add a comma between elements
  //}
  // }
  // std::cout << std::endl;
  
  NumericVector ages(n);
  double step = pfCycleLength / n;
  for (int i = 0; i < n; ++i) {
    ages[i] = (i + 1) * step;  // Assign values directly
  }
  
  NumericVector ys = yfx(ages,inflec);
  
  arma::vec startI0All = beta_starts(betaShape, offset, I0, n);
  
  NumericVector startI0AllNumeric = wrap(startI0All);
  
  NumericVector startI0(2 * n); 
  
  for (int i = 0; i < n; ++i) {
    startI0[i] = ys[i] * startI0AllNumeric[i];
  }
  for (int i = n; i < 2 * n; ++i) {
    startI0[i] = (1 - ys[i - n]) * startI0AllNumeric[i - n];
  } 
  
  std::vector<double> x0(startI0.begin(), startI0.end());
  //   std::cout << "x0: ";
  //   for (size_t i = 0; i < x0.size(); ++i) {
  //      std::cout << x0[i];
  //      if (i < x0.size() - 1) {
  //       std::cout << ", ";  // Add a comma between elements
  //       }
  //  }
  //   std::cout << std::endl;
  
  // Running ODE simulation
  double max_t = 120;
  double dt = 0.1;
  arma::mat odeint_output = constPMR_gammaN_ode_cpp(x0, parmValues, max_t, dt);
  
  double step_size = 30;
  NumericMatrix subsetMatrix = subsetRows(odeint_output, step_size);
  // Rcpp::Rcout << "subsetMatrix: " << subsetMatrix << std::endl;
  
  // Get number of rows
  int numRows = subsetMatrix.nrow();
  NumericVector circ_iRBC_unique(numRows, 0.0); 
  
  for (int j = 1; j <= n; ++j) { 
    circ_iRBC_unique += subsetMatrix(_, j); 
  }
  //Rcpp::Rcout << "circ_iRBC_unique: " << circ_iRBC_unique << std::endl;
  
  
  NumericVector circ_iRBC_rep = repeat_subvector(circ_iRBC_unique);
  // Rcpp::Rcout << "circ_iRBC_rep: " << circ_iRBC_rep << std::endl;
  
  
  double ring_last_stage = round(ring_duration / 24 * n) + 1;
  //Rcpp::Rcout << "ring_last_stage: " << ring_last_stage << std::endl;
  
  if (ring_last_stage <= 2){
    ring_last_stage = 3;
  } else{
    ring_last_stage = ring_last_stage;
  }
  
  NumericVector circ_ring_tot(numRows, 0.0);
  
  for (int j = 1; j < ring_last_stage; ++j) { 
    circ_ring_tot += subsetMatrix(_, j); 
  }
  
  // Rcpp::Rcout << "circ_ring_tot: " << circ_ring_tot << std::endl;
  
  NumericVector ring_prop_estim = circ_ring_tot/circ_iRBC_unique;
  // Rcpp::Rcout << "ring_prop_estim: " << ring_prop_estim << std::endl;
  NumericVector ring_prop_rep = repeat_subvector(ring_prop_estim);
  // Rcpp::Rcout << "ring_prop_rep: " << ring_prop_rep << std::endl;
  
  
  // Calculate SSE
  NumericVector circ_data = data["Circ"]; 
  NumericVector transformed_data(circ_data.size());
  
  for (int i = 0; i < circ_data.size(); ++i) {
    if (NumericVector::is_na(circ_data[i])) {
      transformed_data[i] = NA_REAL;
    } else {
      transformed_data[i] = log10(circ_data[i] + 1.0);
    }
  } 
  
  // Rcpp::Rcout << "transformed_data: " << transformed_data<< std::endl;
  
  NumericVector circ_pred = log10(circ_iRBC_rep + 1.0);
  // Rcpp::Rcout << "circ_pred: " << circ_pred << std::endl;
  
  NumericVector squared_diff_iRBC(transformed_data.size());
  for (int i = 0; i < circ_data.size(); ++i) {
    if (NumericVector::is_na(transformed_data[i])) {
      // If circ_data[i] is NA, set squared_diff[i] to 0
      squared_diff_iRBC[i] = 0.0;
    } else {
      // Otherwise, calculate the squared difference
      squared_diff_iRBC[i] = pow(transformed_data[i] - circ_pred[i], 2);
    }
  } 
  
  double sse_iRBC = sum(squared_diff_iRBC);
  //Rcpp::Rcout << "sse_iRBC: " << sse_iRBC << std::endl;
  
  NumericVector ring_data = data["ring_prop"]; 
  
  //Rcpp::Rcout << "ring_data: " << ring_data << std::endl;
  
  NumericVector squared_diff_ring(ring_data.size());
  for (int i = 0; i < ring_data.size(); ++i) {
    if (NumericVector::is_na(ring_data[i])) {
      squared_diff_ring[i] = 0.0;
    } else {
      squared_diff_ring[i] = pow(ring_data[i] - ring_prop_rep[i], 2);
    }
  } 
  
  double sse_ring = sum(squared_diff_ring);
  //Rcpp::Rcout << "sse_ring: " << sse_ring << std::endl;
  
  double sse = sse_iRBC + sse_ring;
  
  return sse;
  
  
}
