#ifndef DIAGONAL_MPO_INVERSION_HEADER
#define DIAGONAL_MPO_INVERSION_HEADER

#include "../../mps/mps.h"
#include "../../mps/mpo.h"
#include "../../mps/observables.h"
#include "../../dtensor/tensor_index.h"
#include "../../dtensor/tensor_index_op.h"
#include "../../dtensor/dtensor.h"
#include "../../linalg/tensor_cg.h"

template <typename T>
void invdmpo (MPO<T>& H, MPO<T>& invH, int Nsweep, double cutoff, double tol);

template <typename T>
void buildEnv(MPS<T>& psi, MPS<T>& phi, MPO<T>& H, std::vector<dtensor<T>>& MR, std::vector<dtensor<T>>& ML, std::vector<dtensor<T>>& VR, std::vector<dtensor<T>>& VL);

template <typename T>
void updateSite(const char& direction, MPS<T>& psi, MPS<T>& phi, MPO<T>& H, std::vector<dtensor<T>>& MR, std::vector<dtensor<T>>& ML, std::vector<dtensor<T>>& VR, std::vector<dtensor<T>>& VL, const int& site, T& delta);

template <typename T>
void updateEnv(const char& direction, MPS<T>& psi, MPS<T>& phi, MPO<T>& H, std::vector<dtensor<T>>& MR, std::vector<dtensor<T>>& ML, std::vector<dtensor<T>>& VR, std::vector<dtensor<T>>& VL, const int& site);


#endif
