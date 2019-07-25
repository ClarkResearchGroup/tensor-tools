#ifndef DENSE_TENSOR_BASED_DAVIDSON_SOLVER_HEADER
#define DENSE_TENSOR_BASED_DAVIDSON_SOLVER_HEADER

#include "../linalg/lapack_wrapper.h"

#include "../dtensor/dtensor_index.h"
#include "../dtensor/dtensor_index_op.h"
#include "../dtensor/dtensor.h"
//#include "../dtensor/dtensor_view.h"

#include "../qtensor/qtensor_index.h"
#include "../qtensor/qtensor_index_op.h"
#include "../qtensor/qtensor.h"

#include "../qstensor/qstensor.h"

#include "../dtensor/big_dtensor.h"
#include "../qtensor/big_qtensor.h"
#include "../qstensor/big_qstensor.h"

//#include "tensor_cg.h"

template <typename T>
void elemWiseDivide(dtensor<T>& A, double theta, dtensor<T>& B);

template <typename T>
void elemWiseDivide(qtensor<T>& A, double theta, qtensor<T>& B);

template <typename T>
void elemWiseDivide(qstensor<T>& A, double theta, qstensor<T>& B);

// Restarted Davidson method to find the smallest/largest (algebraically) eigenvalue
template <typename T, template <typename> class BigTensorType, template <typename> class TensorType>
double tensor_davidson(BigTensorType<T>& A, TensorType<T>& x, int m, int max_restart, double tol, char mode='S');

#endif
