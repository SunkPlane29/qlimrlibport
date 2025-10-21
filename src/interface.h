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
// - eps_c: Central energy density (in MeV/fm^3).
// - out: Pointer to an output array of length 2 where the results will be stored.
// The output array will contain the mass-radius relation, with the first element being the radius and
// the second element being the mass of the neutron star.
void qlimr_getMR(double *eos_p, double *eos_eps, int length, double eps_c, double R_start, double *out);

// Exported C-compatible function
// This function computes the mass-radius diagram for a neutron star given an equation of state.
// Parameters:
// - eos_p: Pointer to an array of pressure values (in MeV/fm^3).
// - eos_eps: Pointer to an array of energy density values (in MeV/fm^3).
// - length: Length of the eos_p and eos_eps arrays.
// - epsc_start: Starting central energy density (in MeV/fm^3).
// - epsc_end: Ending central energy density (in MeV/fm^3).
// - nstars: Number of stars to compute in the mass-radius diagram.
// - out_epsc: Pointer to an output array where the central energy densities will be stored.
// - out_M: Pointer to an output array where the masses will be stored.
// - out_R: Pointer to an output array where the radii will be stored.
void qlimr_getMRdiagram(double *eos_p, double *eos_eps, int length, double epsc_start, double epsc_end, int nstars, double R_start, double *out_epsc, double *out_M, double *out_R);

#ifdef __cplusplus
}
#endif