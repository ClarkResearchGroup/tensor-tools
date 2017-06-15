#ifndef My_DMRG_DIAG_CLASS_H
#define My_DMRG_DIAG_CLASS_H

#include "../../linalg/arpack_wrapper.h"
#include "../../mps/mps.h"
#include "../../mps/mpo.h"
#include "../../mps/observables.h"
#include "../../dtensor/tensor_index.h"
#include "../../dtensor/tensor_index_op.h"
#include "../../dtensor/dtensor.h"


// "matrix"-"vector" multplication for ARPACK
void arpack_mv(int n, double *in, double *out, MPS<double>& psi, MPO<double>& H, std::vector< dtensor<double> >& TR, std::vector< dtensor<double> >& TL, int site);
void arpack_mv(int n, std::complex<double> *in, std::complex<double> *out, MPS< std::complex<double> >& psi, MPO< std::complex<double> >& H, std::vector< dtensor< std::complex<double> > >& TR, std::vector< dtensor< std::complex<double> > >& TL, int site);

// ARPACK ROUTINES -- Smallest real eigenvalues
void arpack_S(int n, int nev, double *Evals, double **Evecs, double *initVec, MPS<double>& psi, MPO<double>& H, std::vector< dtensor<double> >& TR, std::vector< dtensor<double> >& TL, int site);
void arpack_S(int n, int nev, std::complex<double> *Evals, std::complex<double> **Evecs, std::complex<double> *initVec, MPS< std::complex<double> >& psi, MPO< std::complex<double> >& H, std::vector< dtensor< std::complex<double> > >& TR, std::vector< dtensor< std::complex<double> > >& TL, int site);

// ARPACK ROUTINES -- Largest real eigenvalues
void arpack_L(int n, int nev, double *Evals, double **Evecs, double *initVec, MPS<double>& psi, MPO<double>& H, std::vector< dtensor<double> >& TR, std::vector< dtensor<double> >& TL, int site);
void arpack_L(int n, int nev, std::complex<double> *Evals, std::complex<double> **Evecs, std::complex<double> *initVec, MPS< std::complex<double> >& psi, MPO< std::complex<double> >& H, std::vector< dtensor< std::complex<double> > >& TR, std::vector< dtensor< std::complex<double> > >& TL, int site);

#endif
