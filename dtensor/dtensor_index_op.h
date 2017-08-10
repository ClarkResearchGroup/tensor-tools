#ifndef INDEX_OP_CLASS_FOR_DENSE_TENSOR_HEADER
#define INDEX_OP_CLASS_FOR_DENSE_TENSOR_HEADER

#include "dtensor_index.h"

void find_index_permutation(vector<dtensor_index>& from, vector<dtensor_index>& to, vector<unsigned>& perm);

void index_sets_union(const vector<dtensor_index>& Ac, const vector<dtensor_index>& Bc, vector<dtensor_index>& res);

void index_sets_difference(const vector<dtensor_index>& Ac, const vector<dtensor_index>& Bc, vector<dtensor_index>& res);

void index_sets_intersection(const vector<dtensor_index>& Ac, const vector<dtensor_index>& Bc, vector<dtensor_index>& res);

#endif
