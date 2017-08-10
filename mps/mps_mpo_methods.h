#ifndef MPS_MPO_METHODS_HEADER
#define MPS_MPO_METHODS_HEADER

#include "tt.h"
#include "qtt.h"

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

#endif
