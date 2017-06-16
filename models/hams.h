#ifndef My_HAMS_H
#define My_HAMS_H

#include "../mps/mpo.h"

template <typename T>
void buildHeisenberg(MPO<T>& H, double tE=0);

template <typename T>
void buildHeisenberg(MPO<T>& H, double tE, double* dJ, double* dh);

template <typename T>
void buildSz(MPO<T>& H, int site);

#endif
