/*
 * Copyright 2020 Ryan Levy, Xiongjie Yu, and Bryan K. Clark
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
//recreated ITensor version
template <typename T, template <typename> class BigTensorType, template <typename> class TensorType>
double tensor_davidsonIT(BigTensorType<T>& A, TensorType<T>& x, int m, int max_restart, double tol, char mode='S');

#endif
