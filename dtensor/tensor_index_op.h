#ifndef INDEX_OP_CLASS_FOR_DENSE_TENSOR_HEADER
#define INDEX_OP_CLASS_FOR_DENSE_TENSOR_HEADER

#include "tensor_index.h"

void index_sets_union(const vector<tensor_index>& Ac, const vector<tensor_index>& Bc, vector<tensor_index>& res);

void index_sets_difference(const vector<tensor_index>& Ac, const vector<tensor_index>& Bc, vector<tensor_index>& res);

void index_sets_intersection(const vector<tensor_index>& Ac, const vector<tensor_index>& Bc, vector<tensor_index>& res);

#endif
