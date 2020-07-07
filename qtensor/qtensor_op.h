/*
 * Copyright 2020 Xiongjie Yu, Ryan Levy, and Bryan K. Clark
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *     http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

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

template <typename T>
void svd(qtensor<T>& A,
         vector<qtensor_index>& left, vector<qtensor_index>& right,
         qtensor<T>& U, qtensor<T>& V, vector<double>& S,
         int direction);

template <typename T>
void svd(qtensor<T>& A,
        vector<qtensor_index>& left, vector<qtensor_index>& right,
        qtensor<T>& U, qtensor<T>& V, vector<double>& S,
        int direction, double cutoff);

template <typename T>
void svd(qtensor<T>& A,
        vector<qtensor_index>& left, vector<qtensor_index>& right,
        qtensor<T>& U, qtensor<T>& V, vector<double>& S,
        int direction, double cutoff, long unsigned K);

template <typename T>
void svd_bond(qtensor<T>& A_left, qtensor<T>& A_right,
        qtensor_index& mid, vector<double>& S,
        int direction);

template <typename T>
void svd_bond(qtensor<T>& A_left, qtensor<T>& A_right,
        qtensor_index& mid, vector<double>& S,
        int direction, double cutoff, long unsigned K);

template <typename T>
void svd_bond(qtensor<T>& combined, qtensor<T>& A_left, qtensor<T>& A_right,
        qtensor_index& mid, vector<double>& S,
        int direction, double cutoff, long unsigned K);

#endif
