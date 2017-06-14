#ifndef My_MPO_CLASS
#define My_MPO_CLASS

#include "mpo.h"

template <typename T>
MPO<T>::MPO(int l, int pd, int bd) : TensorTrain<T>() {
	phy_size = pd;
	TensorTrain<T>::setLength(l);
	TensorTrain<T>::setIndexSize(pd*pd);
	TensorTrain<T>::setBondDim(bd);
	TensorTrain<T>::allocateTensors();
}
template MPO<double>::MPO(int l, int pd, int bd);
template MPO< std::complex<double> >::MPO(int l, int pd, int bd);

template <typename T>
MPO<T>::MPO (const MPO<T>& other) : TensorTrain<T>() {
	phy_size = other.phy_size;
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
template MPO<double>::MPO (const MPO<double>& other);
template MPO< std::complex<double> >::MPO (const MPO< std::complex<double> >& other);

template <typename T>
MPO<T>::MPO (MPO<T>&& other) : TensorTrain<T>() {
	phy_size = other.phy_size;
	TensorTrain<T>::setLength(other.length);
	TensorTrain<T>::setIndexSize(other.index_size);
	TensorTrain<T>::setBondDims(other.bond_dims);
	TensorTrain<T>::M = other.M;
	other.M = nullptr;
	TensorTrain<T>::tensors_allocated = true;
	other.tensors_allocated = false;
}
template MPO<double>::MPO (MPO<double>&& other);
template MPO< std::complex<double> >::MPO (MPO< std::complex<double> >&& other);

template <typename T>
void MPO<T>::setMPO(int l, int pd, int bd) {
	phy_size = pd;
	TensorTrain<T>::freeTensors();
	TensorTrain<T>::setLength(l);
	TensorTrain<T>::setIndexSize(pd*pd);
	TensorTrain<T>::setBondDim(bd);
	TensorTrain<T>::allocateTensors();
}
template void MPO<double>::setMPO(int l, int pd, int bd);
template void MPO< std::complex<double> >::setMPO(int l, int pd, int bd);

template <typename T>
MPO<T>::~MPO() {
	TensorTrain<T>::freeTensors();
}
template MPO<double>::~MPO();
template MPO< std::complex<double> >::~MPO();

template <typename T>
MPO<T>& MPO<T>::operator = (const TensorTrain<T>& other){
  assert(other.tensors_allocated);
  if(this!=&other){
		phy_size = int(std::round(std::sqrt(other.index_size)));
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
template MPO<double>& MPO<double>::operator = (const TensorTrain<double>& other);
template MPO< std::complex<double> >& MPO< std::complex<double> >::operator = (const TensorTrain< std::complex<double> >& other);

template <typename T>
void MPO<T>::setIdentity() {
	assert(tensors_allocated);
	TensorTrain<T>::setZero(1);
	for(int i = 0; i < TensorTrain<T>::length; ++i){
		for(int j = 0; j < phy_size; ++j) {
			TensorTrain<T>::M[i][j*phy_size+j].setIdentity();
		}
	}
}
template void MPO<double>::setIdentity();
template void MPO< std::complex<double> >::setIdentity();

template <typename T>
void MPO<T>::square() {
	assert(tensors_allocated);
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> ** tM = new Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> * [TensorTrain<T>::length];
	for(int i = 0; i < TensorTrain<T>::length; ++i){
		tM[i] = new Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> [phy_size*phy_size];
	}
	for(int i = 0; i < TensorTrain<T>::length; ++i){
		for(int j = 0; j < phy_size*phy_size; ++j){
			tM[i][j] = TensorTrain<T>::M[i][j];
		}
	}
	for(int i = 0; i < TensorTrain<T>::length+1; ++i){
		TensorTrain<T>::bond_dims[i] = TensorTrain<T>::bond_dims[i]*TensorTrain<T>::bond_dims[i];
	}
	TensorTrain<T>::setZero();
	for(int i = 0; i < TensorTrain<T>::length; i++){
		int r = tM[i][0].rows();
		int c = tM[i][0].cols();
		for(int p1 = 0; p1 < phy_size; p1++){
			for(int p3 = 0; p3 < phy_size; p3++){
				for(int j = 0; j < r; j++){
					for(int k = 0; k < c; k++){
						for(int p2 = 0; p2 < phy_size; p2++){
							TensorTrain<T>::M[i][p1*phy_size+p3].block(j*r,k*c,r,c) += tM[i][p1*phy_size+p2](j,k) * tM[i][p2*phy_size+p3];
						}
					}
				}
			}
		}
	}
	for(int i = 0; i < TensorTrain<T>::length; ++i){
		delete [] tM[i];
	}
	delete [] tM;
}
template void MPO<double>::square();
template void MPO< std::complex<double> >::square();

template <typename T>
dtensor<T> MPO<T>::tensorize(int site){
	assert(TensorTrain<T>::tensors_allocated);
	assert(site>=0 && site<TensorTrain<T>::length);
	int_vec idx_sizes;
	str_vec idx_names;
	typ_vec idx_types;
	int_vec idx_levels;
	string Link_name_pref = "MPO"+to_string(TensorTrain<T>::TTID)+"Link";
	string Site_name = "Site"+to_string(site);
	if(site==0){
		idx_sizes.push_back(TensorTrain<T>::bond_dims[site+1]);
		idx_sizes.push_back(phy_size);
		idx_sizes.push_back(phy_size);
		idx_names.push_back(Link_name_pref+to_string(site+1));
		idx_names.push_back(Site_name);
		idx_names.push_back(Site_name);
		idx_types.push_back(Link);
		idx_types.push_back(Site);
		idx_types.push_back(Site);
		idx_levels.push_back(0);
		idx_levels.push_back(0);
		idx_levels.push_back(1);
	}else if(site==TensorTrain<T>::length-1){
		idx_sizes.push_back(TensorTrain<T>::bond_dims[site]);
		idx_sizes.push_back(phy_size);
		idx_sizes.push_back(phy_size);
		idx_names.push_back(Link_name_pref+to_string(site));
		idx_names.push_back(Site_name);
		idx_names.push_back(Site_name);
		idx_types.push_back(Link);
		idx_types.push_back(Site);
		idx_types.push_back(Site);
		idx_levels.push_back(0);
		idx_levels.push_back(0);
		idx_levels.push_back(1);
	}else if(site>0 && site<TensorTrain<T>::length-1){
		idx_sizes.push_back(TensorTrain<T>::bond_dims[site]);
		idx_sizes.push_back(TensorTrain<T>::bond_dims[site+1]);
		idx_sizes.push_back(phy_size);
		idx_sizes.push_back(phy_size);
		idx_names.push_back(Link_name_pref+to_string(site));
		idx_names.push_back(Link_name_pref+to_string(site+1));
		idx_names.push_back(Site_name);
		idx_names.push_back(Site_name);
		idx_types.push_back(Link);
		idx_types.push_back(Link);
		idx_types.push_back(Site);
		idx_types.push_back(Site);
		idx_levels.push_back(0);
		idx_levels.push_back(0);
		idx_levels.push_back(0);
		idx_levels.push_back(1);
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
template dtensor<double> MPO<double>::tensorize(int site);
template dtensor< std::complex<double> > MPO< std::complex<double> >::tensorize(int site);


#endif
