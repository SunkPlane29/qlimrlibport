#pragma once

#ifdef __cplusplus
extern "C" {
#endif

// Exported C-compatible function
// This function computes the mass-radius relation for a neutron star given an equation of state.
// Parameters:
// - eos_p: Pointer to an array of pressure values (in MeV/fm^3).
// - eos_eps: Pointer to an array of energy density values (in MeV/fm^3).
// - length: Length of the eos_p and eos_eps arrays.
// - eps_c: Critical energy density (in MeV/fm^3).
// - R_start: Initial radius (in km) for the integration. A common value is 0.0004
// - out: Pointer to an output array where the results will be stored.
// The output array will contain the mass-radius relation, with the first element being the radius and
// the second element being the mass of the neutron star.
void qlimr_getMR(double* eos_p, double* eos_eps, int length, double eps_c, double R_start, double* out);

#ifdef __cplusplus
}
#endif