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

dtensor_index::dtensor_index(unsigned size){
  _size = size;
  _level = 0;
  _name = "i"+to_string(int(1e8*thread_safe_random_double()+1e5*thread_safe_random_double()));
  _type = Link;
  setTag();
}

dtensor_index::dtensor_index(unsigned size, string name){
  _size = size;
  _level = 0;
  _name = name;
  _type = Link;
  setTag();
}

dtensor_index::dtensor_index(unsigned size, string name, index_type type){
  _size = size;
  _level = 0;
  _name = name;
  _type = type;
  setTag();
}

dtensor_index::dtensor_index(unsigned size, string name, index_type type, unsigned level){
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
  return (_tag<A._tag && _level<A._level);
}

inline bool dtensor_index::similar(const dtensor_index& A){
  return (_size==A._size)&&(_name==A._name)&&(_type==A._type)&&((_level==A._level+1)||(_level+1==A._level));
}

inline bool dtensor_index::similar(dtensor_index& A){
  return (_size==A._size)&&(_name==A._name)&&(_type==A._type)&&((_level==A._level+1)||(_level+1==A._level));
}

inline bool dtensor_index::similar(const dtensor_index& A, index_type type){
  return (_size==A._size)&&(_name==A._name)&&(_type==A._type)&&(_type==type)&&((_level==A._level+1)||(_level+1==A._level));
}

inline bool dtensor_index::similar(dtensor_index& A, index_type type){
  return (_size==A._size)&&(_name==A._name)&&(_type==A._type)&&(_type==type)&&((_level==A._level+1)||(_level+1==A._level));
}

#endif
