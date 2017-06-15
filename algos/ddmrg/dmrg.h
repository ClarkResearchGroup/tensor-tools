#ifndef DMRG_CLASS_HEADER
#define DMRG_CLASS_HEADER

#include "dmrg_diag.h"
#include "../../mps/mps.h"
#include "../../mps/mpo.h"
#include "../../mps/observables.h"
#include "../../dtensor/tensor_index.h"
#include "../../dtensor/tensor_index_op.h"
#include "../../dtensor/dtensor.h"

// There are 3 sub routines

template <typename T>
void buildEnv(MPS<T>& psi, MPO<T>& H, std::vector<dtensor<T>>& TR, std::vector<dtensor<T>>& TL);

template <typename T>
void updateSite(const char& direction, MPS<T>& psi, MPO<T>& H, std::vector<dtensor<T>>& TR, std::vector<dtensor<T>>& TL, const int& site, T& EN, char whichEn);

template <typename T>
void updateEnv(const char& direction, MPS<T>& psi, MPO<T>& H, std::vector<dtensor<T>>& TR, std::vector<dtensor<T>>& TL, const int& site);

template <typename T>
T dmrg(MPS<T>& psi, MPO<T>& H, int Nsweep, char whichEn='L');

#endif
