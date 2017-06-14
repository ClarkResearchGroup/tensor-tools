#ifndef My_OBSERVABLES_H
#define My_OBSERVABLES_H

#include "../dtensor/tensor_index.h"
#include "../dtensor/tensor_index_op.h"
#include "../dtensor/dtensor.h"
#include "mps.h"
#include "mpo.h"

//-------------------------------------
// Overlaps
template <typename T>
T psiHphi(MPS<T>& psi, MPO<T>& H, MPS<T>& phi);

template <typename T>
T psiphi(MPS<T>& psi, MPS<T>& phi);


//-------------------------------------
// EE
template <typename T>
void EE(MPS<T>& psi, std::vector<double>& vec);

template <typename T>
void EE(MPO<T>& HH, std::vector<double>& vec);


//-------------------------------------
// MPO operations
template <typename T>
MPO<T> diagonal(const MPO<T>& H);

template <typename T>
MPO<T> offdiagonal(const MPO<T>& H);

template <typename T>
double var(const MPO<T>& H);

template <typename T>
double trace(const MPO<T>& H);

template <typename T>
double l2norm(const MPO<T>& H);

#endif
