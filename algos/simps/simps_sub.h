#ifndef MY_SIMPS_SUB_ROUTINES_H
#define MY_SIMPS_SUB_ROUTINES_H

#include "../../mps/observables.h"
#include "simps_cg.h"

void inv_buildR(MPS& psi, MPO& H, Mxd ** CRM);
void buildMatrix(MPS& psi, MPO& H, Mxd ** CR, Mxd ** CL, const int& site, Mxd& A);
void buildDiag(MPS& psi, MPO& H, Mxd ** CR, Mxd ** CL, const int& site, Mxd& A);
void inv_updateSite(const char& direction, MPS& psi, MPO& H, MPO& HS, Mxd ** CR1, Mxd ** CL1, Mxd ** CR2, Mxd ** CL2, const int& site);
void inv_updateLR(const char& direction, MPS& psi, MPO& H, Mxd ** CR, Mxd ** CL, const int& site);

#endif
