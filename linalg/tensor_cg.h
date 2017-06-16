#ifndef DENSE_TENSOR_BASED_CG_SOLVER_HEADER
#define DENSE_TENSOR_BASED_CG_SOLVER_HEADER

#include "../dtensor/tensor_index.h"
#include "../dtensor/tensor_index_op.h"
#include "../dtensor/dtensor.h"


template <typename T>
void elemWiseDivide(dtensor<T>& A, dtensor<T>& B);

template <typename T>
void LRWx(dtensor<T>& LEnv, dtensor<T>& REnv, dtensor<T>& W, dtensor<T>& x, dtensor<T>& y);

// CG with Jacobi Preconditioner
template <typename T>
int tensor_CG(dtensor<T>& LEnv, dtensor<T>& REnv, dtensor<T>& W, dtensor<T>& x, dtensor<T>& u, int max_iter, double tol, bool preconditioned=false);

#endif
