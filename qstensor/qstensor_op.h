#ifndef QUANTUM_NUMBERED_SPARSE_TENSOR_OPERATIONS_HEADER
#define QUANTUM_NUMBERED_SPARSE_TENSOR_OPERATIONS_HEADER

#include "../util/types_and_headers.h"
#include "../linalg/lapack_wrapper.h"
#include "../qtensor/qtensor_index.h"
#include "../qtensor/qtensor_index_op.h"
#include "qstensor.h"

#define MoveFromLeft 0
#define MoveFromRight 1

template <typename T>
void qr(qstensor<T>& A,
        vector<qtensor_index>& left, vector<qtensor_index>& right,
        qstensor<T>& Q, qstensor<T>& R);

/*template <typename T>
void svd(qstensor<T>& A,
         vector<qtensor_index>& left, vector<qtensor_index>& right,
         qstensor<T>& U, qstensor<T>& V, vector<double>& S,
         int direction);

template <typename T>
void svd(qstensor<T>& A,
        vector<qtensor_index>& left, vector<qtensor_index>& right,
        qstensor<T>& U, qstensor<T>& V, vector<double>& S,
        int direction, double cutoff);*/

template <typename T>
void svd(qstensor<T>& A,
        vector<qtensor_index>& left, vector<qtensor_index>& right,
        qstensor<T>& U, qstensor<T>& V, qtensor<T>& S,
        int direction, double cutoff=0, unsigned K=0);

/*template <typename T>
void svd_bond(qstensor<T>& A_left, qstensor<T>& A_right,
        qtensor_index& mid, vector<double>& S,
        int direction);*/

template <typename T>
void svd_bond(qstensor<T>& A_left, qstensor<T>& A_right,
        qtensor_index& mid, qtensor<T>& S,
        int direction, double cutoff=0, long unsigned K=0);

template <typename T>
void svd_bond(qstensor<T>& combined, qstensor<T>& A_left, qstensor<T>& A_right,
        qtensor_index& mid, qtensor<T>& S,
        int direction, double cutoff=0, long unsigned K=0);

#endif
