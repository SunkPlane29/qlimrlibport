#include "input.h"

#include <iostream>
#include <string>
#include <vector>

Input_QLIMR::Input_QLIMR(double *eos_p, double *eos_eps, int length) {
  params.epsvec.push_back(0.0); // Add zero to the beginning of epsvec
  params.pvec.push_back(0.0);   // Add zero to the beginning of pvec
    for (int i = 0; i < length; ++i) {
        params.pvec.push_back(adimensionalize(eos_p[i], "MeV/fm^3"));
        params.epsvec.push_back(adimensionalize(eos_eps[i], "MeV/fm^3"));
    }
}

// --------------- Adimensionalize conversion function ------------------------
double Input_QLIMR::adimensionalize(double value, std::string unit) {
  double factor;

  if (unit == "MeV/fm^3") {
    factor = 1.0 / 346933.783551;
  } else if (unit == "km") {
    factor = 1.0 / 1.47663;
  } else if (unit == "g/cm^3") {
    factor = 1.0 / 6.17625e+17;
  } else if (unit == "MeV"){
    factor = 8.96162e-61;   
  } else if (unit == "1/fm^3"){
    factor = 3.216297e54;  
  } else if (unit == "-") {
    factor = 1.0;
  } else {
    std::cout << "no unit match to adimensionalize" << std::endl;
    exit(0);
  }
  
  return factor * value;
}
// ----------------------------------------------------------------------------

// --------------- Dimensionalize conversion function ------------------------
double Input_QLIMR::dimensionalize(double value, std::string unit) {
  double factor;

  if (unit == "MeV/fm^3") {
    factor = 346933.783551;
  } else if (unit == "km") {
    factor = 1.47663;
  } else if (unit == "g/cm^3") {
    factor = 6.17625e+17;
  } else if (unit == "MeV"){
    factor = 1.0/8.96162e-61;  
  } else if (unit == "1/fm^3"){
    factor = 1.0 / 3.216297e54;  
  } else if (unit == "Hz") {
    factor = (299792/1.47663); // c [km/s] / lsun [km]
  } else if (unit == "-") { 
    factor = 1.0;
  } else {
    std::cout << "no unit match to dimensionalize" << std::endl;
    exit(0);
  }

  return factor * value;
}
// ----------------------------------------------------------------------------

EOSinterpolation EOS::EoS; // Static structure for EoS

EOS::EOS(gsl_interp_type *type, double *eos_p, double *eos_eps, int length) {
    params = Input_QLIMR(eos_p, eos_eps, length).params;
    // Interpolate p(ε): initialize GSL interpolation spline
    EoS.p_of_e.initialize(type, params.epsvec, params.pvec);

    // std::cout << "epsvec[0] = " << params.epsvec.front() << ", epsvec.back() = " << params.epsvec.back() << "\n";
    // std::cout << "pvec[0] = " << params.pvec.front() << ", pvec.back() = " << params.pvec.back() << "\n";

    // Compute interpolated ε(h) and p(h) 
    calculate_eos_of_h(&params.epsvec, type);
}

// ---------------------- Interpolation class methods -------------------------
void Interpolation::initialize(gsl_interp_type *interp_type, std::vector<double> x,
                               std::vector<double> y) {
  type = interp_type;
  // Size of the independent variable vector                            
  size = x.size();    

  // Allocating GSL interpolation accelerator
  acc = gsl_interp_accel_alloc();

  // Allocating GSL spline of specified type and size
  spline = gsl_spline_alloc(type, size);

  // Initializing GSL spline with given data points (x, y) and size
  gsl_spline_init(spline, x.data(), y.data(), size);
}

// Function to calculate interpolated y value for a given x using GSL spline
double Interpolation::yofx(double x) {
  if (x < 0.0) {
    return 0.0;
  } else if (x > spline->x[spline->size - 1]) {
    return spline->y[spline->size - 1];
  }

  return gsl_spline_eval(spline, x, acc);
}

// Function to calculate dy/dx of interpolated data using GSL spline
double Interpolation::dyofx(double x) {
  if (x < 0.0) {
    return 0.0;
  } else if (x > spline->x[spline->size - 1]) {
    return gsl_spline_eval_deriv(spline, spline->x[spline->size - 1], acc);
  }

  return gsl_spline_eval_deriv(spline, x, acc);
}

// Method to release memory of spline and accelerator
void Interpolation::free() {
  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);
}
//-----------------------------------------------------------------------------

// ---------------------------- Defining dh/dε --------------------------------
double enthalpy_integrand(double epsilon, void *params) {

  double p = ((Interpolation *)params)->yofx(epsilon);
  double cs2 = ((Interpolation *)params)->dyofx(epsilon);

  // std::cout << "epsilon = " << epsilon << ", p = " << p << ", cs2 = " << cs2 << "\n";

  return cs2 / (p + epsilon);
}
// ----------------------------------------------------------------------------

//------------------- Integration: Finding ε(h) and p(h) ----------------------
void EOS::calculate_eos_of_h(std::vector<double> *epsilon,
                                 gsl_interp_type *type) {
  double delta_h;
  double error;
  double h = 0.0;

  // Integration workspace size
  size_t integration_workspace_size = 100 * (epsilon->size());

  // GSL variable type for integrating function
  gsl_function F;

  // Function to be integrated
  F.function = &enthalpy_integrand;

  // Passing object parameter which contains p(ε)
  F.params = &EoS.p_of_e;

  //NOTE: EoS is a static member, so we should be careful with leakages
  EoS.h_vec.clear();
  EoS.e_vec.clear();
  EoS.p_vec.clear();

  // Add h=0 and ε=0 according to definition of h
  EoS.h_vec.push_back(0.0);
  EoS.e_vec.push_back(0.0);
  EoS.p_vec.push_back(0.0);

  // Allocate memory for GSL integration workspace
  gsl_integration_workspace *w =
  gsl_integration_workspace_alloc(integration_workspace_size);

  for (size_t i = 1; i < epsilon->size() - 1; i++) {


    // Perform adaptive quadrature integration using GSL library
    gsl_integration_qag(&F,
                       (*epsilon)[i],     // Lower integration limit
                       (*epsilon)[i + 1], // Upper integration limit
                       1e-9,              // Absolute tolerance
                       1e-9,              // Relative tolerance
                       1000,              // Maximal number of subintervals
                       6,                 // Integration method (Gauss-Kronrod)
                       w,                 // Work space
                       &delta_h,          // Estimated step size
                       &error             // Estimated integration error
    );

    // std::cout << h + delta_h << "    " << delta_h << "\n";

    //NOTE: When delta_h is too small, there will be an interpolation error because two h's will have the same value
    if (delta_h < 5e-7) {
      // std::cout << "Zero delta_h encountered at i = " << i << ", epsilon = " << (*epsilon)[i] << "\n";
      break;
    }

    h = h + delta_h;
    EoS.h_vec.push_back(h);
    EoS.e_vec.push_back((*epsilon)[i + 1]);
    EoS.p_vec.push_back(EoS.p_of_e.yofx((*epsilon)[i + 1]));
  }

  gsl_integration_workspace_free(w);

  EoS.h_of_e.initialize(type, EoS.e_vec, EoS.h_vec); // Interpolate h(ε)
  EoS.h_of_p.initialize(type, EoS.p_vec, EoS.h_vec); // Interpolate h(p)
  EoS.e_of_h.initialize(type, EoS.h_vec, EoS.e_vec); // Interpolate ε(h)
  EoS.p_of_h.initialize(type, EoS.h_vec, EoS.p_vec); // Interpolate p(h)
}
//-----------------------------------------------------------------------------

EOS::~EOS() {
  EoS.e_of_h.free();
  EoS.p_of_h.free();
  EoS.h_of_e.free();
  EoS.h_of_p.free();
  EoS.p_of_e.free();
}