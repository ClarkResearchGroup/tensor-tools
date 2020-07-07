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


#ifndef DENSE_TENSOR_OPERATIONS_HEADER
#define DENSE_TENSOR_OPERATIONS_HEADER

#include "../util/types_and_headers.h"
#include "../linalg/lapack_wrapper.h"
#include "dtensor_index.h"
#include "dtensor_index_op.h"
#include "dtensor.h"
//#include "dtensor_view.h"

string indToStr(vector<dtensor_index> &indices,unordered_map<string,char> &charMap);
string indicesToChar(vector<dtensor_index> &indices, unordered_map<string,char> &charMap);
unsigned indicesToSize(vector<dtensor_index> &indices);

template <typename T>
void qr(dtensor<T>& A,
        vector<dtensor_index>& left, vector<dtensor_index>& right,
        dtensor<T>& Q, dtensor<T>& R);

/*template <typename T>
void qr(dtensor_view<T>& A,
        vector<dtensor_index>& left, vector<dtensor_index>& right,
        dtensor<T>& Q, dtensor<T>& R);*/

/*template <typename T>
void svd(dtensor<T>& A,
         vector<dtensor_index>& left, vector<dtensor_index>& right,
         dtensor<T>& U, dtensor<T>& V, vector<double>& S,
         int direction);

  template <typename T>
void svd(dtensor<T>& A,
    vector<dtensor_index>& left,
    vector<dtensor_index>& right,
    dtensor<T>& U, dtensor<T>& V, dtensor<T>& S,
    int direction);*/

/*template <typename T>
void svd(dtensor_view<T>& A,
        vector<dtensor_index>& left, vector<dtensor_index>& right,
        dtensor<T>& U, dtensor<T>& V, vector<double>& S,
        int direction);*/

/*template <typename T>
void svd(dtensor<T>& A,
        vector<dtensor_index>& left, vector<dtensor_index>& right,
        dtensor<T>& U, dtensor<T>& V, vector<double>& S,
        int direction, double cutoff);*/

/*template <typename T>
void svd(dtensor_view<T>& A,
        vector<dtensor_index>& left, vector<dtensor_index>& right,
        dtensor<T>& U, dtensor<T>& V, vector<double>& S,
        int direction, double cutoff);*/

template <typename T>
void svd(dtensor<T>& A,
        vector<dtensor_index>& left, vector<dtensor_index>& right,
        dtensor<T>& U, dtensor<T>& V, vector<double>& S,
        int direction, double cutoff=0, long unsigned K=0);

/*template <typename T>
void svd(dtensor_view<T>& A,
        vector<dtensor_index>& left, vector<dtensor_index>& right,
        dtensor<T>& U, dtensor<T>& V, vector<double>& S,
        int direction, double cutoff, long unsigned K);*/

template <typename T>
void svd(dtensor<T>& A,
    vector<dtensor_index>& left,
    vector<dtensor_index>& right,
    dtensor<T>& U, dtensor<T>& V, dtensor<T>& S,
    int direction,double cutoff=0, long unsigned K=0);

/*template <typename T>
void svd_bond(dtensor<T>& A_left, dtensor<T>& A_right,
        dtensor_index& mid, vector<double>& S,
        int direction);*/

/*template <typename T>
void svd_bond(dtensor_view<T>& A_left, dtensor_view<T>& A_right,
        dtensor_index& mid, vector<double>& S,
        int direction);*/

template <typename T>
void svd_bond(dtensor<T>& A_left, dtensor<T>& A_right,
        dtensor_index& mid, vector<double>& S,
        int direction, double cutoff=0, long unsigned K=0);

  template <typename T>
void svd_bond(dtensor<T>& combined, dtensor<T>& A_left, dtensor<T>& A_right,
    dtensor_index& mid, dtensor<T>& S,
    int direction, double cutoff=0, long unsigned K=0);

/*template <typename T>
void svd_bond(dtensor_view<T>& A_left, dtensor_view<T>& A_right,
        dtensor_index& mid, vector<double>& S,
        int direction, double cutoff, long unsigned K);*/

template <typename T>
void svd_bond(dtensor<T>& combined, dtensor<T>& A_left, dtensor<T>& A_right,
        dtensor_index& mid, vector<double>& S,
        int direction, double cutoff=0, long unsigned K=0);

/*template <typename T>
void svd_bond(dtensor<T>& combined, dtensor_view<T>& A_left, dtensor_view<T>& A_right,
        dtensor_index& mid, vector<double>& S,
        int direction, double cutoff, long unsigned K);*/

#endif
