#ifndef ALL_MPS_RELATED_HEADERS
#define ALL_MPS_RELATED_HEADERS

#include "tt.h"
#include "qtt.h"
#include "qstt.h"
#include "observables.h"
#include "mps_mpo_methods.h"
namespace ezh5{
	template<> Node& Node::operator = (QN_t_1 val){
		if(id==-1){
			hid_t dataspace_id = H5Screate(H5S_SCALAR);
			this->id = H5Dcreate(pid, path.c_str(), H5T_NATIVE_INT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		}
		hid_t error_id = H5Dwrite(id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &val.qn);
		H5Dclose(this->id);
		this->id = -1;
		return *this;
	}
	template<> Node& Node::operator >> (QN_t_1& val){
		hid_t dataset_id  = H5Dopen(pid, this->path.c_str(), H5P_DEFAULT);
		hid_t datatype_id = H5Dget_type(dataset_id);
		hid_t error_id = H5Dread(dataset_id, datatype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, &val.qn);
		H5Dclose(dataset_id);
		return *this;
	}
	template<> Node& Node::operator = (std::vector<QN_t_1>& vec){
		hid_t dataspace_id = -1;
		if(this->id == -1){
			hsize_t dims[1];
			dims[0] = vec.size();
			dataspace_id = H5Screate_simple(1, dims, NULL);
			assert(dataspace_id>=0);
			this->id = H5Dcreate(pid, path.c_str(), H5T_NATIVE_INT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		}
    vector<int> temp;temp.reserve(vec.size()); for(auto& c: vec) temp.push_back(c.qn);
		hid_t error_id = H5Dwrite(id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &temp[0]);
		H5Dclose(this->id);
		this->id = -1;
		if (dataspace_id != -1) {H5Sclose(dataspace_id);}
		return *this;
	};
	template<> Node& Node::operator >> (std::vector<QN_t_1>& vec){
		hid_t dataset_id = H5Dopen2(pid, this->path.c_str(), H5P_DEFAULT);
		hid_t datatype_id = H5Dget_type(dataset_id);
		hid_t dataspace_id = H5Dget_space(dataset_id);
		hsize_t dims[1];
		int err = H5Sget_simple_extent_dims(dataspace_id,dims,NULL); assert(err>=0);
    vector<int> temp;
		if (dims[0]>0) {
				temp.resize(dims[0]);
				hid_t error_id = H5Dread(dataset_id, datatype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, &temp[0]);
				H5Dclose(dataset_id);
		}
    for(auto c: temp) vec.push_back({c});
		return *this;
	};
}

#endif
