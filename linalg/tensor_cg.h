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

#ifndef DENSE_TENSOR_BASED_CG_SOLVER_HEADER
#define DENSE_TENSOR_BASED_CG_SOLVER_HEADER

#include "../dtensor/dtensor_index.h"
#include "../dtensor/dtensor_index_op.h"
#include "../dtensor/dtensor.h"
#include "../dtensor/dtensor_view.h"
#include "../dtensor/big_dtensor.h"

#include "../qtensor/qtensor_index.h"
#include "../qtensor/qtensor_index_op.h"
#include "../qtensor/qtensor.h"
#include "../qtensor/big_qtensor.h"

template <typename T>
void elemWiseDivide(dtensor<T>& A, dtensor<T>& B);

template <typename T>
void elemWiseDivide(qtensor<T>& A, qtensor<T>& B);

// CG, optional Jacobi Preconditioner
template <typename T, template <typename> class BigTensorType, template <typename> class TensorType1, template <typename> class TensorType2>
int tensor_CG(BigTensorType<T>& A, TensorType1<T>& x, TensorType2<T>& b, int max_iter, double tol, bool preconditioning=false);

#endif
