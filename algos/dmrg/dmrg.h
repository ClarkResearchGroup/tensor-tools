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
