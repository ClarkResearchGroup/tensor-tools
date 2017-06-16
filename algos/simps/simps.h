#ifndef SHIFT_AND_INVERSE_MPS_HEADER
#define SHIFT_AND_INVERSE_MPS_HEADER

#include "../../mps/mps.h"
#include "../../mps/mpo.h"
#include "../../mps/observables.h"
#include "../../dtensor/tensor_index.h"
#include "../../dtensor/tensor_index_op.h"
#include "../../dtensor/dtensor.h"
#include "../../linalg/tensor_cg.h"

template <typename T>
void simps (double tE, MPO<T>& H, MPS<T>& psi, int Nsweep);

template <typename T>
void buildEnv(MPS<T>& psi, MPO<T>& H, MPO<T>& HSq, std::vector<dtensor<T>>& MR, std::vector<dtensor<T>>& ML, std::vector<dtensor<T>>& VR, std::vector<dtensor<T>>& VL);

template <typename T>
void updateSite(const char& direction, MPS<T>& psi, MPO<T>& H, MPO<T>& HSq, std::vector<dtensor<T>>& MR, std::vector<dtensor<T>>& ML, std::vector<dtensor<T>>& VR, std::vector<dtensor<T>>& VL, const int& site, T& En, T& EnSq);

template <typename T>
void updateEnv(const char& direction, MPS<T>& psi, MPO<T>& H, MPO<T>& HSq, std::vector<dtensor<T>>& MR, std::vector<dtensor<T>>& ML, std::vector<dtensor<T>>& VR, std::vector<dtensor<T>>& VL, const int& site);


#endif
