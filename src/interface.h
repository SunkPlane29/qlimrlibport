#pragma once

#ifdef __cplusplus
extern "C" {
#endif

// Exported C-compatible function
void qlimr_getMR(double* eos_p, double* eos_eps, int length, double eps_c, double R_start, double* out);

#ifdef __cplusplus
}
#endif