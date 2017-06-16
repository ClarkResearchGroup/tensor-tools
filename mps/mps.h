#ifndef My_MPS_CLASS_H
#define My_MPS_CLASS_H

#include "tt.h"
#include "../dtensor/tensor_index.h"
#include "../dtensor/tensor_index_op.h"
#include "../dtensor/dtensor.h"

/*
Most of the MPS methods are inherited from the TensorTrain base class, like
(1) Storage management
(2) Canonicalization, SVD;
(3) Basic arithmetic operations (+,-,*,/)
*/

template <typename T>
class MPS : public TensorTrain<T>
{
public:
	MPS () : TensorTrain<T>() {};
	MPS (int l, int pd, int bd);
	MPS (const MPS<T>& other);
	MPS (MPS<T>&& other);
	~ MPS ();

	void setMPS(int l, int pd, int bd);
	void checkBonds();

	double maxBondDim();
	double avgBondDim();

	MPS& operator = (const TensorTrain<T>& other);
	MPS& operator = (const MPS<T>& other);

	dtensor<T> tensorize(int site);
};

#endif
