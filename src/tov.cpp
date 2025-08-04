#include "tov.h"
#include "input.h"

#include <cmath>
#include <iostream>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

int TOV_equations(double h, const double y[], double f[], void *paraM_sol)
{

  // Recast input parameter pointer of the ODE system as an input object
  EOSinterpolation *eos = (EOSinterpolation *)paraM_sol;

  // EoS variables ε(h) and p(h) to be inserted into the ODE system
  double e = eos->e_of_h.yofx(h);
  double p = eos->p_of_h.yofx(h);

  // Functions R(h) and M(h) to be solved
  double R = y[0];
  double M = y[1];

  double x0 = 4 * M_PI * std::pow(R, 3);

  // TOV equations in h-formulation form. Here, f[0]= dR/dh and f[1]= dM/dh.
  f[0] = -R * (-2 * M + R) / (M + x0 * p);
  f[1] = -e * x0 * (-2 * M + R) / (M + p * x0);

  return GSL_SUCCESS;
}

// Initial conditions for TOVh system 
Initial_conditions_TOV IC_TOV(double epsilon_c, EOSinterpolation &EoS, double R_start)
{

  // Values at exactly the center
  double pc = EoS.p_of_e.yofx(epsilon_c);
  double hc = EoS.h_of_e.yofx(epsilon_c);
  double Cs2 = EoS.p_of_e.dyofx(epsilon_c);

//   std::cout << "pc = " << pc << ", hc = " << hc << ", Cs2 = " << Cs2 << "\n"; 

  // ################## Asymptotic solutions ###################

  double R = R_start; //NOTE: hardcoded

  double M = (4.0 / 3.0) * M_PI * epsilon_c * pow(R, 3);

  // -----------------------------------------------------------

  double p = pc - (2.0 / 3.0) * M_PI * pow(R, 2) * (epsilon_c + 3.0 * pc) *
                      (epsilon_c + pc);

  double e = epsilon_c - (2.0 / 3.0) * M_PI * pow(R, 2) *
                             (epsilon_c + 3.0 * pc) * (epsilon_c + pc) / Cs2;

  // ####################### dM/dh, dR/dh at R=Rε ########################

  double dRdh = -(R * (R - 2.0 * M)) / (M + 4.0 * M_PI * pow(R, 3) * p);

  double dMdh = 4.0 * M_PI * e * pow(R, 2) * dRdh;

  auto IC_tov = Initial_conditions_TOV();
  // Filling initial conditions and initial step size into structure
  IC_tov.R_start = R;
  IC_tov.M_start = M;
  IC_tov.h_start = hc - (2.0 / 3.0) * M_PI * pow(R, 2) * (epsilon_c + 3.0 * p);
  IC_tov.h_istep = 0.1 * std::min(abs(R / dRdh), abs(M / dMdh));

  // cout << min(abs(R / dRdh), abs(M / dMdh)) << endl;
//   std::cout << "hc" << hc << std::endl;
  // //exit(1);

  return IC_tov;
}
