#ifndef INVDMPO_CG_SOLVER_H
#define INVDMPO_CG_SOLVER_H

#include "../../mps/observables.h"

Mxd simps_mult(MPO& H, Mxd& v, Mxd ** CRM, Mxd ** CLM, const int& site);

// CG with Jacobi Preconditioner -- (efficient for Diagonal dominant matrix)
int CG_linearSolver(int n, MPO& H, Mxd& M, Mxd& v, Mxd& b, Mxd ** CRM, Mxd ** CLM, const int& site);

#endif
