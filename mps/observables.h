#ifndef My_OBSERVABLES_H
#define My_OBSERVABLES_H

#include "tt.h"
#include "qtt.h"
#include "qstt.h"

//-------------------------------------
// Overlaps
template <typename T>
T psiHphi(MPS<T>& psi, MPO<T>& H, MPS<T>& phi);

template <typename T>
T psiHphi(qMPS<T>& psi, qMPO<T>& H, qMPS<T>& phi);

template <typename T>
T psiHphi(qsMPS<T>& psi, qsMPO<T>& H, qsMPS<T>& phi);

template <typename T>
T psiHKphi(MPS<T>& psi, MPO<T>& H, MPO<T>& K, MPS<T>& phi);

template <typename T>
T psiHKphi(qMPS<T>& psi, qMPO<T>& H, qMPO<T>& K, qMPS<T>& phi);

template <typename T>
T psiphi(MPS<T>& psi, MPS<T>& phi);

template <typename T>
T psiphi(qMPS<T>& psi, qMPS<T>& phi);

template <typename T>
T psiphi(qsMPS<T>& psi, qsMPS<T>& phi);

#endif
