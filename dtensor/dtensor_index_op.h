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

#ifndef INDEX_OP_CLASS_FOR_DENSE_TENSOR_HEADER
#define INDEX_OP_CLASS_FOR_DENSE_TENSOR_HEADER

#include "dtensor_index.h"

void find_index_permutation(vector<dtensor_index>& from, vector<dtensor_index>& to, vector<unsigned>& perm);

void index_sets_union(const vector<dtensor_index>& Ac, const vector<dtensor_index>& Bc, vector<dtensor_index>& res);

void index_sets_difference(const vector<dtensor_index>& Ac, const vector<dtensor_index>& Bc, vector<dtensor_index>& res);

void index_sets_intersection(const vector<dtensor_index>& Ac, const vector<dtensor_index>& Bc, vector<dtensor_index>& res);

#endif
