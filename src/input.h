#pragma once

#include <string>
#include <vector>

#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>

struct Input_params {
    std::vector<double> pvec;
    std::vector<double> epsvec;
};

class Input_QLIMR {
public:
    Input_params params;

    Input_QLIMR() = default;
    Input_QLIMR(double *eos_p, double *eos_eps, int length);

    // Unit conversion functions
    static double adimensionalize(double value, std::string unit);
    static double dimensionalize(double value, std::string unit);
};

class Interpolation : public Input_QLIMR {
public:
    // GSL-type variables for interpolation
    gsl_interp_accel *acc;
    gsl_spline       *spline;
    gsl_interp_type  *type;

    // Size of data to be interpolated
    size_t size;

    // Initializing GSL spline for interpolation
    void initialize(gsl_interp_type *type, std::vector<double> x, std::vector<double> y);

    // Evaluate interpolated function at point x
    double yofx(double x);

    // Evaluate derivative of interpolated function at point x
    double dyofx(double x);

    // Free GSL spline and accelerator memory
    void free();
};

struct EOSinterpolation {

    // Interpolation type objects:

    Interpolation p_of_e; // Pressure as a function of energy density
    Interpolation h_of_e; // Pseudo-enthalpy as a function of energy density
    Interpolation h_of_p; // Enthalpy as a function of pressure
    Interpolation e_of_h; // Energy density as a function of pseudo-enthalpy
    Interpolation p_of_h; // Pressure as a function of pseudo-enthalpy

    // Vectors to store EoS data columns

    std::vector<double> e_vec;  // Energy density vector
    std::vector<double> p_vec;  // Pressure vector
    std::vector<double> h_vec;  // Enthalpy vector
};

class EOS : public Interpolation {
public:
  // Static structure for EoS to be used along the code 
  static EOSinterpolation EoS;

  EOS() = default;
  // Parametric constructor: 
  EOS(gsl_interp_type *type, double *eos_p, double *eos_eps, int length);

  ~EOS();

  // Method function to compute EoS in terms of pseudo-enthalpy (h)
  void calculate_eos_of_h(std::vector<double> *epsilon, gsl_interp_type *type);
};