#ifndef DENSE_TENSOR_BASED_CG_SOLVER_HEADER
#define DENSE_TENSOR_BASED_CG_SOLVER_HEADER

#include "../dtensor/dtensor_index.h"
#include "../dtensor/dtensor_index_op.h"
#include "../dtensor/dtensor.h"
#include "../dtensor/dtensor_view.h"

#include "../qtensor/qtensor_index.h"
#include "../qtensor/qtensor_index_op.h"
#include "../qtensor/qtensor.h"

template <typename T>
void elemWiseDivide(dtensor<T>& A, dtensor<T>& B);

template <typename T>
void elemWiseDivide(qtensor<T>& A, qtensor<T>& B);

// CG, optional Jacobi Preconditioner
template <typename T, template <typename> class BigTensorType, template <typename> class TensorType1, template <typename> class TensorType2>
int tensor_CG(BigTensorType<T>& A, TensorType1<T>& x, TensorType2<T>& b, int max_iter, double tol, bool preconditioning=false);

#endif
