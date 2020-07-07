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


#ifndef EZH5
#define EZH5

#include "ezh5.h"

namespace ezh5{
	//Define Data Types
	template<typename T> struct DataType{
		static hid_t id;
	};

	// complex number datatype
	template<typename T> struct DataType<std::complex<T> >{
		static hid_t id;
		DataType(){
	    	H5Tinsert(id, "r", 0, DataType<T>::id);  // the name 'r' is to be compatible with h5py
	    	H5Tinsert(id, "i", sizeof(T), DataType<T>::id);
		}
  };

  template<> hid_t DataType<bool>::id = H5T_NATIVE_HBOOL;
  template<> hid_t DataType<char>::id = H5T_NATIVE_CHAR;
  template<> hid_t DataType<double>::id = H5T_NATIVE_DOUBLE;
  template<> hid_t DataType<float>::id = H5T_NATIVE_FLOAT;
  template<> hid_t DataType<int>::id = H5T_NATIVE_INT;
  template<> hid_t DataType<long>::id = H5T_NATIVE_LONG;
  template<> hid_t DataType<unsigned int>::id = H5T_NATIVE_UINT;
  template<> hid_t DataType<unsigned long>::id = H5T_NATIVE_ULONG;
  template<> hid_t DataType<std::complex<int> >::id = H5Tcreate(H5T_COMPOUND, sizeof(std::complex<int>));
  template<> hid_t DataType<std::complex<float> >::id = H5Tcreate(H5T_COMPOUND, sizeof(std::complex<float>));
  template<> hid_t DataType<std::complex<double> >::id = H5Tcreate(H5T_COMPOUND, sizeof(std::complex<double>));

  DataType<std::complex<float> > obj_to_run_constructor_float;
  DataType<std::complex<double> > obj_to_run_constructor_double;
}

namespace ezh5{
	///////////////////////////
	Node Node::operator[] (const std::string& path_more){
		if(this->id==-1){
			htri_t is_exist = H5Lexists(pid, path.c_str(), H5P_DEFAULT);
			if (is_exist<0){
				assert(false);
			}else if (is_exist==false){
				this->id = H5Gcreate2(pid, path.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
			}else{
				this->id = H5Gopen(pid, path.c_str(), H5P_DEFAULT);
			}
		}
		assert(this->id>0);
		return Node(this->id, path_more, verbose_);
	}
	////////
	template<typename T> Node& Node::operator = (T val){
		if(id==-1){
			hid_t dataspace_id = H5Screate(H5S_SCALAR);
			this->id = H5Dcreate(pid, path.c_str(), DataType<T>::id, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		}
		hid_t error_id = H5Dwrite(id, DataType<T>::id, H5S_ALL, H5S_ALL, H5P_DEFAULT, &val);
		H5Dclose(this->id);
		this->id = -1;
		return *this;
	}
	template Node& Node::operator = (bool val);
	template Node& Node::operator = (int val);
  template Node& Node::operator = (unsigned val);
  template Node& Node::operator = (unsigned long val);
  template Node& Node::operator = (long val);
	template Node& Node::operator = (float val);
	template Node& Node::operator = (double val);
	template Node& Node::operator = (char val);
	template Node& Node::operator = (std::complex<int> val);
	template Node& Node::operator = (std::complex<float> val);
	template Node& Node::operator = (std::complex<double> val);
	////////
	template<typename T> Node& Node::operator >> (T& val){
		hid_t dataset_id  = H5Dopen(pid, this->path.c_str(), H5P_DEFAULT);
		hid_t datatype_id = H5Dget_type(dataset_id);
		hid_t error_id = H5Dread(dataset_id, datatype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, &val);
		H5Dclose(dataset_id);
		return *this;
	}
	template Node& Node::operator >> (bool& val);
	template Node& Node::operator >> (int& val);
	template Node& Node::operator >> (long& val);
  template Node& Node::operator >> (unsigned& val);
  template Node& Node::operator >> (unsigned long& val);
	template Node& Node::operator >> (float& val);
	template Node& Node::operator >> (double& val);
	template Node& Node::operator >> (char& val);
	template Node& Node::operator >> (std::complex<int>& val);
	template Node& Node::operator >> (std::complex<float>& val);
	template Node& Node::operator >> (std::complex<double>& val);
	////////
	template<typename T> Node& Node::operator = (std::vector<T>& vec){
		hid_t dataspace_id = -1;
		if(this->id == -1){
			hsize_t dims[1];
			dims[0] = vec.size();
			dataspace_id = H5Screate_simple(1, dims, NULL);
			assert(dataspace_id>=0);
			this->id = H5Dcreate(pid, path.c_str(), DataType<T>::id, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		}
		hid_t error_id = H5Dwrite(id, DataType<T>::id, H5S_ALL, H5S_ALL, H5P_DEFAULT, &vec[0]);
		H5Dclose(this->id);
		this->id = -1;
		if (dataspace_id != -1) {H5Sclose(dataspace_id);}
		return *this;
	}
	template Node& Node::operator = (std::vector<int>& vec);
	template Node& Node::operator = (std::vector<long>& vec);
  template Node& Node::operator = (std::vector<unsigned>& vec);
  template Node& Node::operator = (std::vector<unsigned long>& vec);
	template Node& Node::operator = (std::vector<float>& vec);
	template Node& Node::operator = (std::vector<double>& vec);
	template Node& Node::operator = (std::vector<char>& vec);
	template Node& Node::operator = (std::vector<std::complex<int> >& vec);
	template Node& Node::operator = (std::vector<std::complex<float> >& vec);
	template Node& Node::operator = (std::vector<std::complex<double> >& vec);
	////////
	template<typename T> Node& Node::operator >> (std::vector<T>& vec){
		hid_t dataset_id = H5Dopen2(pid, this->path.c_str(), H5P_DEFAULT);
		hid_t datatype_id = H5Dget_type(dataset_id);
		hid_t dataspace_id = H5Dget_space(dataset_id);
		hsize_t dims[1];
		int err = H5Sget_simple_extent_dims(dataspace_id,dims,NULL); assert(err>=0);
		if (dims[0]>0) {
				vec.resize(dims[0]);
				hid_t error_id = H5Dread(dataset_id, datatype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, &vec[0]);
				H5Dclose(dataset_id);
		}
		return *this;
	}
	template Node& Node::operator >> (std::vector<int>& vec);
	template Node& Node::operator >> (std::vector<long>& vec);
  	template Node& Node::operator >> (std::vector<unsigned>& vec);
  	template Node& Node::operator >> (std::vector<unsigned long>& vec);
	template Node& Node::operator >> (std::vector<float>& vec);
	template Node& Node::operator >> (std::vector<double>& vec);
	template Node& Node::operator >> (std::vector<char>& vec);
	template Node& Node::operator >> (std::vector<std::complex<int> >& vec);
	template Node& Node::operator >> (std::vector<std::complex<float> >& vec);
	template Node& Node::operator >> (std::vector<std::complex<double> >& vec);
	////////
	template<typename T> Node& Node::operator = (Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& m){
		hsize_t dims[2];
		dims[0] = m.rows();
		dims[1] = m.cols();
		if(id==-1){
			hid_t dataspace_id = H5Screate_simple(2, dims, NULL);
			assert(dataspace_id>=0);
			this->id = H5Dcreate(pid, path.c_str(), DataType<T>::id, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
			assert(this->id>=0);
		}
		hid_t error_id = H5Dwrite(id, DataType<T>::id, H5S_ALL, H5S_ALL, H5P_DEFAULT, m.data());
		H5Dclose(this->id);
		this->id = -1;
		return *this;
	}
	template Node& Node::operator = (Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>& m);
  	template Node& Node::operator = (Eigen::Matrix<unsigned, Eigen::Dynamic, Eigen::Dynamic>& m);
	template Node& Node::operator = (Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>& m);
	template Node& Node::operator = (Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& m);
	template Node& Node::operator = (Eigen::Matrix<std::complex<int>, Eigen::Dynamic, Eigen::Dynamic>& m);
	template Node& Node::operator = (Eigen::Matrix<std::complex<float>, Eigen::Dynamic, Eigen::Dynamic>& m);
	template Node& Node::operator = (Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>& m);
	////////
	template<typename T> Node& Node::operator >> (Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& m){
		hid_t dataset_id   = H5Dopen(pid, this->path.c_str(), H5P_DEFAULT);
		hid_t datatype_id  = H5Dget_type(dataset_id);
		hid_t dataspace_id = H5Dget_space(dataset_id);
		hsize_t rank       = H5Sget_simple_extent_ndims(dataspace_id);
		assert(rank==2);
		hsize_t dims[2];
		H5Sget_simple_extent_dims(dataspace_id, dims, NULL);
		m.setZero(dims[0],dims[1]);
		hid_t error_id = H5Dread(dataset_id, datatype_id, H5S_ALL, dataspace_id, H5P_DEFAULT, m.data());
		H5Dclose(dataset_id);
		return *this;
	}
	template Node& Node::operator >> (Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>& m);
  	template Node& Node::operator >> (Eigen::Matrix<unsigned, Eigen::Dynamic, Eigen::Dynamic>& m);
	template Node& Node::operator >> (Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>& m);
	template Node& Node::operator >> (Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& m);
	template Node& Node::operator >> (Eigen::Matrix<std::complex<int>, Eigen::Dynamic, Eigen::Dynamic>& m);
	template Node& Node::operator >> (Eigen::Matrix<std::complex<float>, Eigen::Dynamic, Eigen::Dynamic>& m);
	template Node& Node::operator >> (Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>& m);
	///////////////////////////
	File::File(const std::string& path, unsigned flags, bool verbose){
		verbose_ = verbose;
		if(flags==H5F_ACC_RDWR || flags==H5F_ACC_RDONLY){
			this->id = H5Fopen(path.c_str(), flags, H5P_DEFAULT);
		}else if (flags==H5F_ACC_TRUNC || flags==H5F_ACC_EXCL){
			this->id = H5Fcreate(path.c_str(), flags, H5P_DEFAULT, H5P_DEFAULT);
		}else{
			assert(false && "unknow file access mode");
		}
	}
	///////////////////////////
}


#endif
