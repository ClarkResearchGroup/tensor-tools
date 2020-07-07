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

#ifndef MPS_MPO_METHODS_HEADER
#define MPS_MPO_METHODS_HEADER

#include "tt.h"
#include "qtt.h"
#include "observables.h"

template <typename T, unsigned N>
void sum(vector< dTensorTrain<T, N> >& x, dTensorTrain<T, N>& res, int max_iter, double cutoff, int max_bd, double tol=1e-12, bool verbose=false);

template <typename T, unsigned N>
void sum(vector< qTensorTrain<T, N> >& x, qTensorTrain<T, N>& res, int max_iter, double cutoff, int max_bd, double tol=1e-12, bool verbose=false);

template <typename T>
void mult(MPO<T>& A, MPO<T>& B, MPO<T>& res, int max_iter, double cutoff, int max_bd, double tol=1e-12, bool verbose=false);

template <typename T>
void mult(qMPO<T>& A, qMPO<T>& B, qMPO<T>& res, int max_iter, double cutoff, int max_bd, double tol=1e-12, bool verbose=false);

template <typename T>
void mult(MPO<T>& A, MPS<T>& B, MPS<T>& res, int max_iter, double cutoff, int max_bd, double tol=1e-12, bool verbose=false);

template <typename T>
void mult(qMPO<T>& A, qMPS<T>& B, qMPS<T>& res, int max_iter, double cutoff, int max_bd, double tol=1e-12, bool verbose=false);

template<typename T>
MPS<T> exactApplyMPO(MPO<T> & K, MPS<T> & psi,double cutoff=1E-13,int maxm=-1, bool verbose=false);

template <typename T>  
T errorMPOProd(MPS<T> & psi2, MPO<T> & K, MPS<T> & psi1);

#endif
