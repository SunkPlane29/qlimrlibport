#include "interface.h"
#include "input.h"
#include "tov.h"

#include <iostream>
#include <cmath>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

extern "C" void qlimr_getMR(double *eos_p, double *eos_eps, int length, double eps_c, double *out) {

    // std::cout << "p[0] = " << eos_p[0] << ", p[end] = " << eos_p[length - 1] << "\n";
    // std::cout << "eps[0] = " << eos_eps[0] << ", eps[end] = " << eos_eps[length - 1] << "\n";
    // Number of dependent variables of the ODE system to be solved
    const int dim = 2;

    // std::cout << "getting eos" << "\n";
    EOS eos((gsl_interp_type *)gsl_interp_steffen, eos_p, eos_eps, length);

    // print first and last of h vec
    // std::cout << "h_vec[0] = " << eos.EoS.h_vec.front() << ", h_vec.back() = " << eos.EoS.h_vec.back() << "\n";
    // std::cout << "e_vec[0] = " << eos.EoS.e_vec.front() << ", e_vec.back() = " << eos.EoS.e_vec.back() << "\n";
    // std::cout << "p_vec[0] = " << eos.EoS.p_vec.front() << ", p_vec.back() = " << eos.EoS.p_vec.back() << "\n";

    // std::cout << "eos obtained" << "\n";

    // Defining GSL variables: system, step, control and evolve 
    gsl_odeiv2_system sys = {TOV_equations, NULL, dim, &eos.EoS};
    gsl_odeiv2_step *s = gsl_odeiv2_step_alloc(gsl_odeiv2_step_rkf45, dim);
    gsl_odeiv2_control *c = gsl_odeiv2_control_y_new(1e-10, 1e-10);
    gsl_odeiv2_evolve *e = gsl_odeiv2_evolve_alloc(dim);

    // std::cout << "gsl variables allocated" << "\n";

    std::cout << "eps_c = " << eps_c << "\n"; 
    eps_c = EOS::adimensionalize(eps_c, "MeV/fm^3");
    // std::cout << "eps_c = " << eps_c << "\n";

    // Obtain initial conditions
    Initial_conditions_TOV IC_tov = IC_TOV(eps_c, eos.EoS);

    // std::cout << "Initial conditions obtained: R_start = " << IC_tov.R_start
            //   << ", M_start = " << IC_tov.M_start
            //   << ", h_start = " << IC_tov.h_start
            //   << ", h_istep = " << IC_tov.h_istep << "\n";

    // Initial conditions array for r and m at R_start
    double y[2] = {IC_tov.R_start, IC_tov.M_start};

    // Initial value of h to begin the integration
    double h = IC_tov.h_start;

    // Initial step size for the GSL adaptative step size algorithm
    double h_istep = -1.0 * IC_tov.h_istep;

    // Final value of h to stop integration
    double h_stop = 0.0;

    // Solution vectors to store h, R, M, p, e and ν at each step in h
    std::vector<double> h_sol, R_sol, M_sol, p_sol, e_sol, nu_sol;

    // Add initial conditions at h_start
    h_sol.push_back(h);
    R_sol.push_back(y[0]);
    M_sol.push_back(y[1]);
    p_sol.push_back(eos.EoS.p_of_e.yofx((eos.EoS.e_of_h.yofx(h))));
    e_sol.push_back(eos.EoS.e_of_h.yofx(h));

    // std::cout << "Initial conditions set: R = " << y[0] << ", M = " << y[1] << ", h = " << h << "\n";

    // Integration loop from h = h_start up to the surface when h = 0.0 
    while (h > h_stop) {

        // std::cout << "Integrating at h = " << h << "\n";

        // Solve the ODE system using GSL with rkf45 method at each step
        int status = gsl_odeiv2_evolve_apply(e, c, s,
            &sys, &h, h_stop, &h_istep, y);

        // std::cout << "After integration: R = " << y[0] << ", M = " << y[1] << ", h = " << h << "\n";

        // std::cout << "After integration: R = " << y[0] << ", M = " << y[1] << ", h = " << h << "\n";

        // Stop solving if something's wrong with the numerical integration
        if (status != GSL_SUCCESS) break;
            
        // Store h values to be used for finding ν(h) */
        h_sol.push_back(h);

        // Store radial distance (R) after each step
        R_sol.push_back(y[0]);

        // Store enclosed mass (M) after each step
        M_sol.push_back(y[1]);

        // Store pressure (p) after each step
        p_sol.push_back(eos.EoS.p_of_e.yofx((eos.EoS.e_of_h.yofx(h))));

        // Store energy density (ε) after each step
        e_sol.push_back(eos.EoS.e_of_h.yofx(h));
    }

    // // Check that integration indeed goes up to h = 0 where p = 0. 
    // std::cout << "h = " << h  << " " << "p = " << p_sol.back() << " " << "e = " << e_sol.back()  << "\n"; 

    // Free GSL variables of the integrator
    gsl_odeiv2_evolve_free(e);
    gsl_odeiv2_control_free(c);
    gsl_odeiv2_step_free(s);

    // Storing total mass M = m(h=0) and radius R = r(h=0) 
    out[0] = M_sol.back();
    out[1] = Input_QLIMR::dimensionalize(R_sol.back(), "km");

    std::cout << "Final mass M = " << out[0] << ", Final radius R = " << out[1] << "\n";
}