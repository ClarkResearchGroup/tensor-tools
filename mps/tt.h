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

#ifndef DENSE_TENSORTRAIN_CLASS_H
#define DENSE_TENSORTRAIN_CLASS_H

#include "../util/ezh5.h"
#include "../linalg/lapack_wrapper.h"
#include "../dtensor/dtensor_index.h"
#include "../dtensor/dtensor_index_op.h"
#include "../dtensor/dtensor.h"
#include "../dtensor/dtensor_view.h"
#include "../dtensor/dtensor_op.h"
#include "../dtensor/big_dtensor.h"
#include "../models/sites/sites.h"

/*
This class is the base class of MPS and MPO.
*/

// N is the number of physical indices per site (1:MPS, 2:MPO)
template <typename T, unsigned N>
class dTensorTrain
{
public:
  //---------------------------------------------------------------------------
  // ID used for naming the the tensors
  unsigned _id; // assigned randomly when an object is created

  //---------------------------------------------------------------------------
  // Basic properties
  unsigned length;         // length of the tensor train
  unsigned phy_dim;        // dimension of each physical bond
  int center;              // canonicalization center
  uint_vec bond_dims;       // dimensions of virtual bonds

  //---------------------------------------------------------------------------
  // Data containers
  vector< dtensor<T> > A;

  //---------------------------------------------------------------------------
  // Flag -- whether the tensors are allocated
  bool tensors_allocated;

  //---------------------------------------------------------------------------
  // Constructors
  dTensorTrain();
  dTensorTrain(abstract_sites* s, unsigned bd = 1);
  dTensorTrain(abstract_sites* s, str_vec product_string);
  dTensorTrain(unsigned l, unsigned pd, unsigned bd = 1);
  dTensorTrain(const dTensorTrain<T, N>& other);
  dTensorTrain(dTensorTrain<T, N>&& other);
  // Destructor
  ~dTensorTrain(){};

  //---------------------------------------------------------------------------
  // Set shapes
  void setLength(int L);
  void setPhysicalDim(int s);
  void setBondDim(int bd);
  void setBondDims(const uint_vec& bds);

  //---------------------------------------------------------------------------
  // Tensor data management
  void allocateTensors(unsigned* product_state=nullptr);
  void freeTensors();

  //---------------------------------------------------------------------------
  // Set tensors values
  void setZero();
  void setRandom();

  //---------------------------------------------------------------------------
  // Basic arithmetic operations
  dTensorTrain<T, N>& operator = (const dTensorTrain<T, N>& other);
  dTensorTrain<T, N>& operator = (dTensorTrain<T, N>&& other);
  dTensorTrain<T, N>& operator *=(const T c);
  dTensorTrain<T, N>& operator /=(const T c);
  dTensorTrain<T, N>  operator * (const T c) const;
  dTensorTrain<T, N>  operator / (const T c) const;

  //---------------------------------------------------------------------------
  // Print dTensorTrain information
  void print(int level=0);

  //---------------------------------------------------------------------------
  // HDF5 storage
  void save(std::string fn);
  void load(std::string fn);
  void save(ezh5::Node& fW);
  void load(ezh5::Node& fR);

  //---------------------------------------------------------------------------
  // Canonicalization
  void rc();
  void lc();

  //---------------------------------------------------------------------------
  // Adjust canonical center
  double position(int site);

  //---------------------------------------------------------------------------
  // Norm
  void normalize();
  double norm();
};


// Definition of MPS and MPO through dTensorTrain
template <typename T>
using MPS =  dTensorTrain<T, 1>;

template <typename T>
using MPO =  dTensorTrain<T, 2>;


#endif
