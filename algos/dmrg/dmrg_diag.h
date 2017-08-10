#ifndef My_DMRG_DIAG_CLASS_H
#define My_DMRG_DIAG_CLASS_H

#include "../../linalg/arpack_wrapper.h"
#include "../../mps/observables.h"

void dmrg_mtv(int n, double *in, double *out, MPS& psi, MPO& H, Mxd ** CRM, Mxd ** CLM, int site);

void dsaupd_LE(int n, int nev, double *Evals, double **Evecs, double *initVec, MPS& psi, MPO& H, Mxd ** CRM, Mxd ** CLM, int site);

void dsaupd_HE(int n, int nev, double *Evals, double **Evecs, double *initVec, MPS& psi, MPO& H, Mxd ** CRM, Mxd ** CLM, int site);

#endif
