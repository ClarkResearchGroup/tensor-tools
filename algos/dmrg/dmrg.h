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

#ifndef DMRG_CLASS_HEADER
#define DMRG_CLASS_HEADER

#include "../../mps/tt.h"
#include "../../mps/qtt.h"
#include "../../mps/observables.h"
#include "../../linalg/tensor_davidson.h"


//-------------------------------------------------------------------------------------------
// MPS, MPO
template <typename T>
void buildEnv(MPS<T>& psi, MPO<T>& H, std::vector<dtensor<T>>& TR, std::vector<dtensor<T>>& TL);
template <typename T>
void updateSite(MPS<T>& psi, MPO<T>& H, std::vector<dtensor<T>>& TR, std::vector<dtensor<T>>& TL, const unsigned& site, T& energy, int& direction, int max_bd, double cutoff, char mode, int search_space_size, int max_restart);
template <typename T>
void updateEnv(MPS<T>& psi, MPO<T>& H, std::vector<dtensor<T>>& TR, std::vector<dtensor<T>>& TL, const unsigned& site, int& direction);
template <typename T>
T dmrg(MPS<T>& psi, MPO<T>& H, int num_sweeps, int max_bd = 100, double cutoff=1e-12, char mode='S', int search_space_size=3, int max_restart=1,int start_sweep=0);

template <typename T>
T dmrg(MPS<T>& psi, MPO<T>& H, int num_sweeps, const std::vector<int>& max_bd, const std::vector<double>& cutoff, const std::vector<int>& max_restart);
//-------------------------------------------------------------------------------------------
// qMPS, qMPO
template <typename T>
void buildEnv(qMPS<T>& psi, qMPO<T>& H, std::vector<qtensor<T>>& TR, std::vector<qtensor<T>>& TL);
template <typename T>
void updateSite(qMPS<T>& psi, qMPO<T>& H, std::vector<qtensor<T>>& TR, std::vector<qtensor<T>>& TL, const unsigned& site, T& energy, int& direction, int max_bd, double cutoff, char mode, int search_space_size, int max_restart);
template <typename T>
void updateEnv(qMPS<T>& psi, qMPO<T>& H, std::vector<qtensor<T>>& TR, std::vector<qtensor<T>>& TL, const unsigned& site, int& direction);
template <typename T>
T dmrg(qMPS<T>& psi, qMPO<T>& H, int num_sweeps, int max_bd = 100, double cutoff=1e-12, char mode='S', int search_space_size=3, int max_restart=1,int start_sweep=0);

template <typename T>
T dmrg(qMPS<T>& psi, qMPO<T>& H, int num_sweeps, const std::vector<int>& max_bd, const std::vector<double>& cutoff, const std::vector<int>& max_restart);
#endif
