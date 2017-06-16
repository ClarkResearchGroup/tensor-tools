#ifndef MPS_MPO_OPERATIONS_HEADER
#define MPS_MPO_OPERATIONS_HEADER

#include "../dtensor/tensor_index.h"
#include "../dtensor/tensor_index_op.h"
#include "../dtensor/dtensor.h"
#include "mps.h"
#include "mpo.h"

// Kronecker product of Eigen matrices
template <typename T>
void KroneckerProd(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& A, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& B, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& res);

// Apply (MPO) A (transposed if needed) to (MPS) psi, send (and compress if needed) the result to (MPS) res
template <typename T>
void fitApplyMPO(MPO<T>& A, MPS<T>& phi, MPS<T>& res, char transA = 'I', bool to_compress=true, double cutoff=1e-10);

// Apply (MPO) A (transposed if needed) to (MPO) B, send (and compress if needed) the result to (MPO) res
template <typename T>
void fitApplyMPO(MPO<T>& A, MPO<T>& B, MPO<T>& res, char transA = 'I', bool to_compress=true, double cutoff=1e-10);

// diagonal of MPO
template <typename T>
MPO<T> diagonal(const MPO<T>& H);

// offdiagonal of MPO
template <typename T>
MPO<T> offdiagonal(const MPO<T>& H);

// take the diagonal of a MPO and pretend it is a MPS
template <typename T>
MPS<T> DiagonalMPOAsMPS(const MPO<T>& H);

// take the diagonal of a MPO and pretend it is a MPS
template <typename T>
MPO<T> MPSAsDiagonalMPO(const MPS<T>& H);

// average variance of MPO
template <typename T>
double var(const MPO<T>& H);

// trace of MPO
template <typename T>
T trace(const MPO<T>& H);

// l2norm of MPO
template <typename T>
double l2norm(const MPO<T>& H);

#endif
