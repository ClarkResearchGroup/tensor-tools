#ifndef QUANTUM_NUMBERED_TENSOR_OPERATIONS_HEADER
#define QUANTUM_NUMBERED_TENSOR_OPERATIONS_HEADER

#include "../util/types_and_headers.h"
#include "../linalg/lapack_wrapper.h"
#include "qtensor_index.h"
#include "qtensor_index_op.h"
#include "qtensor.h"

#define MoveFromLeft 0
#define MoveFromRight 1

template <typename T>
void qr(qtensor<T>& A,
        vector<qtensor_index>& left, vector<qtensor_index>& right,
        qtensor<T>& Q, qtensor<T>& R);

/*template <typename T>
void svd(qtensor<T>& A,
         vector<qtensor_index>& left, vector<qtensor_index>& right,
         qtensor<T>& U, qtensor<T>& V, vector<double>& S,
         int direction);

template <typename T>
void svd(qtensor<T>& A,
        vector<qtensor_index>& left, vector<qtensor_index>& right,
        qtensor<T>& U, qtensor<T>& V, vector<double>& S,
        int direction, double cutoff);*/

template <typename T>
void svd(qtensor<T>& A,
        vector<qtensor_index>& left, vector<qtensor_index>& right,
        qtensor<T>& U, qtensor<T>& V, qtensor<T>& S,
        int direction, double cutoff=0, unsigned K=0);

/*template <typename T>
void svd_bond(qtensor<T>& A_left, qtensor<T>& A_right,
        qtensor_index& mid, vector<double>& S,
        int direction);*/

template <typename T>
void svd_bond(qtensor<T>& A_left, qtensor<T>& A_right,
        qtensor_index& mid, qtensor<T>& S,
        int direction, double cutoff=0, long unsigned K=0);

template <typename T>
void svd_bond(qtensor<T>& combined, qtensor<T>& A_left, qtensor<T>& A_right,
        qtensor_index& mid, qtensor<T>& S,
        int direction, double cutoff=0, long unsigned K=0);

#endif
