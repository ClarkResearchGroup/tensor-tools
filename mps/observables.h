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

#endif
