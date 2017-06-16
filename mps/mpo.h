#ifndef My_MPO_CLASS_H
#define My_MPO_CLASS_H

#include "tt.h"
#include "../dtensor/tensor_index.h"
#include "../dtensor/tensor_index_op.h"
#include "../dtensor/dtensor.h"

/*
Most of the MPO methods are inherited from the TensorTrain base class, like
(1) Storage management
(2) Canonicalization, SVD;
(3) Basic arithmetic operations (+,-,*,/)
*/

template <typename T>
class MPO : public TensorTrain<T>
{
public:
	MPO () : TensorTrain<T>() {};
	MPO (int l, int pd, int bd);
	MPO (const MPO<T>& other);
	MPO (MPO<T>&& other);
	~ MPO ();

	void setMPO(int l, int pd, int bd);
	void checkBonds();

	void setIdentity(); // defaults to bond dimension 1
	void square();

	MPO& operator = (const TensorTrain<T>& other);
	MPO& operator = (const MPO<T>& other);

	dtensor<T> tensorize(int site);
};

#endif
