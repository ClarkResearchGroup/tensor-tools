#ifndef INDEX_OP_CLASS_FOR_QUANTUM_NUMBERED_TENSOR_HEADER
#define INDEX_OP_CLASS_FOR_QUANTUM_NUMBERED_TENSOR_HEADER

#include "qtensor_index.h"

void find_index_permutation(vector<qtensor_index>& from, vector<qtensor_index>& to, vector<unsigned>& perm);

void index_sets_union(const vector<qtensor_index>& Ac, const vector<qtensor_index>& Bc, vector<qtensor_index>& res);

void index_sets_difference(const vector<qtensor_index>& Ac, const vector<qtensor_index>& Bc, vector<qtensor_index>& res);

void index_sets_intersection(const vector<qtensor_index>& Ac, const vector<qtensor_index>& Bc, vector<qtensor_index>& res);

#endif
