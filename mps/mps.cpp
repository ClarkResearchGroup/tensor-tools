#ifndef My_MPS_CLASS
#define My_MPS_CLASS

#include "mps.h"

template <typename T>
MPS<T>::MPS(int l, int pd, int bd) : TensorTrain<T>() {
	TensorTrain<T>::setLength(l);
	TensorTrain<T>::setIndexSize(pd);
	TensorTrain<T>::setBondDim(bd);
	checkBonds();
	TensorTrain<T>::allocateTensors();
	TensorTrain<T>::setRandom();
	TensorTrain<T>::normalize();
}
template MPS<double>::MPS(int l, int pd, int bd);
template MPS< std::complex<double> >::MPS(int l, int pd, int bd);

template <typename T>
MPS<T>::MPS (const MPS<T>& other) : TensorTrain<T>() {
	TensorTrain<T>::setLength(other.length);
	TensorTrain<T>::setIndexSize(other.index_size);
	TensorTrain<T>::setBondDims(other.bond_dims);
	checkBonds();
	TensorTrain<T>::allocateTensors();
	for (int i = 0; i < TensorTrain<T>::length; i++){
		for (int j = 0; j < TensorTrain<T>::index_size; j++) {
			TensorTrain<T>::M[i][j] = other.M[i][j];
		}
	}
}
template MPS<double>::MPS (const MPS<double>& other);
template MPS< std::complex<double> >::MPS (const MPS< std::complex<double> >& other);

template <typename T>
MPS<T>::MPS (MPS<T>&& other) : TensorTrain<T>() {
	TensorTrain<T>::setLength(other.length);
	TensorTrain<T>::setIndexSize(other.index_size);
	TensorTrain<T>::setBondDims(other.bond_dims);
	TensorTrain<T>::M = other.M;
	other.M = nullptr;
	TensorTrain<T>::tensors_allocated = true;
	other.tensors_allocated = false;
}
template MPS<double>::MPS (MPS<double>&& other);
template MPS< std::complex<double> >::MPS (MPS< std::complex<double> >&& other);

template <typename T>
void MPS<T>::setMPS(int l, int pd, int bd) {
	TensorTrain<T>::freeTensors();
	TensorTrain<T>::setLength(l);
	TensorTrain<T>::setIndexSize(pd);
	TensorTrain<T>::setBondDim(bd);
	checkBonds();
	TensorTrain<T>::allocateTensors();
	TensorTrain<T>::setRandom();
	TensorTrain<T>::normalize();
}
template void MPS<double>::setMPS(int l, int pd, int bd);
template void MPS< std::complex<double> >::setMPS(int l, int pd, int bd);

template <typename T>
void MPS<T>::checkBonds(){
	for (int i = 1; i < TensorTrain<T>::length; i++) {
		TensorTrain<T>::bond_dims[i] = std::min(TensorTrain<T>::bond_dims[i], TensorTrain<T>::index_size*TensorTrain<T>::bond_dims[i-1]);
	}
	for (int i = TensorTrain<T>::length-1; i > 0; i--) {
		TensorTrain<T>::bond_dims[i] = std::min(TensorTrain<T>::bond_dims[i], TensorTrain<T>::index_size*TensorTrain<T>::bond_dims[i+1]);
	}
}
template void MPS<double>::checkBonds();
template void MPS< std::complex<double> >::checkBonds();

template <typename T>
MPS<T>::~MPS() {
	TensorTrain<T>::freeTensors();
}
template MPS<double>::~MPS();
template MPS< std::complex<double> >::~MPS();

template <typename T>
MPS<T>& MPS<T>::operator = (const TensorTrain<T>& other){
  assert(other.tensors_allocated);
  if(this!=&other){
		TensorTrain<T>::freeTensors();
    TensorTrain<T>::setLength(other.length);
    TensorTrain<T>::setIndexSize(other.index_size);
    TensorTrain<T>::setBondDims(other.bond_dims);
    TensorTrain<T>::allocateTensors();
    for (int i = 0; i < TensorTrain<T>::length; i++){
      for (int j = 0; j < TensorTrain<T>::index_size; j++) {
        TensorTrain<T>::M[i][j] = other.M[i][j];
      }
    }
  }
  return *this;
}
template MPS<double>& MPS<double>::operator = (const TensorTrain<double>& other);
template MPS< std::complex<double> >& MPS< std::complex<double> >::operator = (const TensorTrain< std::complex<double> >& other);

template <typename T>
double MPS<T>::maxBondDim(){
	assert(TensorTrain<T>::tensors_allocated);
	return *std::max_element(TensorTrain<T>::bond_dims.begin(), TensorTrain<T>::bond_dims.end());
}
template double MPS<double>::maxBondDim();
template double MPS< std::complex<double> >::maxBondDim();

template <typename T>
double MPS<T>::avgBondDim(){
	assert(TensorTrain<T>::tensors_allocated);
	double avgBD = 0.0;
	for(auto v : TensorTrain<T>::bond_dims) avgBD += v;
	return avgBD/(TensorTrain<T>::length+1);
}
template double MPS<double>::avgBondDim();
template double MPS< std::complex<double> >::avgBondDim();

template <typename T>
dtensor<T> MPS<T>::tensorize(int site){
	assert(TensorTrain<T>::tensors_allocated);
	assert(site>=0 && site<TensorTrain<T>::length);
	int_vec idx_sizes;
	str_vec idx_names;
	typ_vec idx_types;
	int_vec idx_levels;
	string Link_name_pref = "MPS"+to_string(TensorTrain<T>::TTID)+"Link";
	string Site_name = "Site"+to_string(site);
	if(site==0){
		idx_sizes.push_back(TensorTrain<T>::bond_dims[site+1]);
		idx_sizes.push_back(TensorTrain<T>::index_size);
		idx_names.push_back(Link_name_pref+to_string(site+1));
		idx_names.push_back(Site_name);
		idx_types.push_back(Link);
		idx_types.push_back(Site);
		idx_levels.push_back(0);
		idx_levels.push_back(0);
	}else if(site==TensorTrain<T>::length-1){
		idx_sizes.push_back(TensorTrain<T>::bond_dims[site]);
		idx_sizes.push_back(TensorTrain<T>::index_size);
		idx_names.push_back(Link_name_pref+to_string(site));
		idx_names.push_back(Site_name);
		idx_types.push_back(Link);
		idx_types.push_back(Site);
		idx_levels.push_back(0);
		idx_levels.push_back(0);
	}else if(site>0 && site<TensorTrain<T>::length-1){
		idx_sizes.push_back(TensorTrain<T>::bond_dims[site]);
		idx_sizes.push_back(TensorTrain<T>::bond_dims[site+1]);
		idx_sizes.push_back(TensorTrain<T>::index_size);
		idx_names.push_back(Link_name_pref+to_string(site));
		idx_names.push_back(Link_name_pref+to_string(site+1));
		idx_names.push_back(Site_name);
		idx_types.push_back(Link);
		idx_types.push_back(Link);
		idx_types.push_back(Site);
		idx_levels.push_back(0);
		idx_levels.push_back(0);
		idx_levels.push_back(0);
	}
	dtensor<T> A(idx_sizes,idx_names,idx_types,idx_levels);
	int k = 0;
	for (int i = 0; i < TensorTrain<T>::index_size; i++) {
		for (int j = 0; j < TensorTrain<T>::M[site][i].size(); j++) {
			A._T.data()[k] = TensorTrain<T>::M[site][i].data()[j];
			++k;
		}
	}
	return A;
}
template dtensor<double> MPS<double>::tensorize(int site);
template dtensor< std::complex<double> > MPS< std::complex<double> >::tensorize(int site);

#endif
