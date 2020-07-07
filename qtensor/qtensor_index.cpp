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


#ifndef INDEX_CLASS_FOR_QUANTUM_NUMBERED_TENSOR
#define INDEX_CLASS_FOR_QUANTUM_NUMBERED_TENSOR

#include "qtensor_index.h"

qtensor_index::qtensor_index(){
  _arrow = Inward;
  _level = 0;
  _name = "qi"+to_string(int(1e8*thread_safe_random_double()+1e6*thread_safe_random_double()));
  _type = Link;
  setTag();
  _sorted=false;
}

qtensor_index::qtensor_index(arrow_type arrow){
  _arrow = arrow;
  _level = 0;
  _name = "qi"+to_string(int(1e8*thread_safe_random_double()+1e5*thread_safe_random_double()));
  _type = Link;
  setTag();
  _sorted=false;
}

qtensor_index::qtensor_index(arrow_type arrow, string name){
  _arrow = arrow;
  _level = 0;
  _name = name;
  _type = Link;
  setTag();
  _sorted=false;
}

qtensor_index::qtensor_index(arrow_type arrow, string name, index_type type){
  _arrow = arrow;
  _level = 0;
  _name = name;
  _type = type;
  setTag();
  _sorted=false;
}

qtensor_index::qtensor_index(arrow_type arrow, string name, index_type type, unsigned level){
  _arrow = arrow;
  _level = level;
  _name = name;
  _type = type;
  setTag();
  _sorted=false;
}

qtensor_index::qtensor_index(const qtensor_index& other){
  _arrow = other._arrow;
  _level = other._level;
  _name = other._name;
  _type = other._type;
  _qn = other._qn;
  _tag = other._tag;
  _sorted = other._sorted;
}

qtensor_index::qtensor_index(qtensor_index&& other){
  _arrow = other._arrow;
  _level = other._level;
  _name = other._name;
  _type = other._type;
  _qn = std::move(other._qn);
  _tag = other._tag;
  _sorted = other._sorted;
}

qtensor_index qtensor_index::operator = (const qtensor_index& other){
  if(this!=&other){
    _arrow = other._arrow;
    _level = other._level;
    _name = other._name;
    _type = other._type;
    _qn = other._qn;
    _tag = other._tag;
    _sorted = other._sorted;
  }
  return *this;
}

bool qtensor_index::operator == (const qtensor_index& A) const{
  return (_tag==A._tag && _arrow==A._arrow && _level==A._level);
}

bool qtensor_index::operator != (const qtensor_index& A) const{
  return (_tag!=A._tag || _arrow!=A._arrow || _level!=A._level);
}

bool qtensor_index::operator <  (const qtensor_index& A) const{
  return (_tag<A._tag && _arrow<A._arrow && _level<A._level);
  //return std::tie(_tag,_arrow,_level) < std::tie(A._tag,A._arrow,A._level); //strict weak ordering
}

qtensor_index noPrime(qtensor_index& index) {
  qtensor_index newIndex = index;
  newIndex.noPrime();
  newIndex.setTag();
  return newIndex;
}
#endif
