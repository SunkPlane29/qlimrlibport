#pragma once

#include "input.h"

struct Initial_conditions_TOV {
    double R_start;
    double M_start; 
    double h_start;
    double h_istep;
  };

int TOV_equations(double h, const double y[], double f[], void *paraM_sol);
Initial_conditions_TOV IC_TOV(double epsilon_c, EOSinterpolation &EoS);