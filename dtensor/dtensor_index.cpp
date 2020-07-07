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

#include "dtensor_index.h"

dtensor_index::dtensor_index(){
  _size = 1;
  _level = 0;
  _name = "i"+to_string(int(1e8*thread_safe_random_double()+1e6*thread_safe_random_double()));
  _type = Link;
  setTag();
}

dtensor_index::dtensor_index(int64_t size){
  _size = size;
  _level = 0;
  _name = "i"+to_string(int(1e8*thread_safe_random_double()+1e5*thread_safe_random_double()));
  _type = Link;
  setTag();
}

dtensor_index::dtensor_index(int64_t size, string name){
  _size = size;
  _level = 0;
  _name = name;
  _type = Link;
  setTag();
}

dtensor_index::dtensor_index(int64_t size, string name, index_type type){
  _size = size;
  _level = 0;
  _name = name;
  _type = type;
  setTag();
}

dtensor_index::dtensor_index(int64_t size, string name, index_type type, unsigned level){
  _size = size;
  _level = level;
  _name = name;
  _type = type;
  setTag();
}

dtensor_index::dtensor_index(const dtensor_index& other){
  _size = other._size;
  _level = other._level;
  _name = other._name;
  _type = other._type;
  _tag = other._tag;
}

dtensor_index dtensor_index::operator = (const dtensor_index& other){
  if(this!=&other){
    _size = other._size;
    _level = other._level;
    _name = other._name;
    _type = other._type;
    _tag = other._tag;
  }
  return *this;
}

bool dtensor_index::operator == (const dtensor_index& A) const{
  // return (_size==A._size)&&(_name==A._name)&&(_type==A._type)&&(_level==A._level);
  return (_tag==A._tag && _level==A._level);
}

bool dtensor_index::operator != (const dtensor_index& A) const{
  // return (_size!=A._size)||(_name!=A._name)||(_type!=A._type)||(_level!=A._level);
  return (_tag!=A._tag || _level!=A._level);
}

bool dtensor_index::operator <  (const dtensor_index& A) const{
  //return (_tag<A._tag && _level<A._level);
  return std::tie(_tag,_level) < std::tie(A._tag,A._level);
}


dtensor_index noPrime(dtensor_index& index) {
  dtensor_index newIndex = index;
  newIndex.noPrime();
  newIndex.setTag();
  return newIndex;
}
#endif
