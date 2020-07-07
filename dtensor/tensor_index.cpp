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


#ifndef INDEX_CLASS_FOR_DENSE_TENSOR
#define INDEX_CLASS_FOR_DENSE_TENSOR

#include "tensor_index.h"

tensor_index::tensor_index(unsigned size){
  _size = size;
  _level = 0;
  _name = "i"+to_string(int(1e8*drand48()+1e5*drand48()));
  _type = Link;
  _tag = to_string(_size)+" "+_name+" "+to_string(_type)+" "+to_string(_level);
}

tensor_index::tensor_index(unsigned size, string name){
  _size = size;
  _level = 0;
  _name = name;
  _type = Link;
  _tag = to_string(_size)+" "+_name+" "+to_string(_type)+" "+to_string(_level);
}

tensor_index::tensor_index(unsigned size, string name, index_type type){
  _size = size;
  _level = 0;
  _name = name;
  _type = type;
  _tag = to_string(_size)+" "+_name+" "+to_string(_type)+" "+to_string(_level);
}

tensor_index::tensor_index(unsigned size, string name, index_type type, unsigned level){
  _size = size;
  _level = level;
  _name = name;
  _type = type;
  _tag = to_string(_size)+" "+_name+" "+to_string(_type)+" "+to_string(_level);
}

tensor_index::tensor_index(const tensor_index& other){
  _size = other._size;
  _level = other._level;
  _name = other._name;
  _type = other._type;
  _tag = other._tag;
}

bool tensor_index::operator == (const tensor_index& A) const{
  // return (_size==A._size)&&(_name==A._name)&&(_type==A._type)&&(_level==A._level);
  return (_tag==A._tag);
}

bool tensor_index::operator != (const tensor_index& A) const{
  // return (_size!=A._size)||(_name!=A._name)||(_type!=A._type)||(_level!=A._level);
  return (_tag!=A._tag);
}

bool tensor_index::operator <  (const tensor_index& A) const{
  return (_tag<A._tag);
}

inline bool tensor_index::similar(const tensor_index& A){
  return (_size==A._size)&&(_name==A._name)&&(_type==A._type)&&(_level!=A._level);
}

inline bool tensor_index::similar(tensor_index& A){
  return (_size==A._size)&&(_name==A._name)&&(_type==A._type)&&(_level!=A._level);
}

inline bool tensor_index::similar(const tensor_index& A, index_type type){
  return (_size==A._size)&&(_name==A._name)&&(_type==A._type)&&(_type==type)&&(_level!=A._level);
}

inline bool tensor_index::similar(tensor_index& A, index_type type){
  return (_size==A._size)&&(_name==A._name)&&(_type==A._type)&&(_type==type)&&(_level!=A._level);
}

#endif
