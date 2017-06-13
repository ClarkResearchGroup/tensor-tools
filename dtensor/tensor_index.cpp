#ifndef INDEX_CLASS_FOR_DENSE_TENSOR
#define INDEX_CLASS_FOR_DENSE_TENSOR

#include "tensor_index.h"

tensor_index::tensor_index(unsigned size){
  _size = size;
  _level = 0;
  _name = "i"+to_string(int(1e8*drand48()+1e5*drand48()));
  _type = Link;
}

tensor_index::tensor_index(unsigned size, string name){
  _size = size;
  _level = 0;
  _name = name;
  _type = Link;
}

tensor_index::tensor_index(unsigned size, string name, index_type type){
  _size = size;
  _level = 0;
  _name = name;
  _type = type;
}

tensor_index::tensor_index(unsigned size, string name, index_type type, unsigned level){
  _size = size;
  _level = level;
  _name = name;
  _type = type;
}

tensor_index::tensor_index(const tensor_index& other){
  _size = other._size;
  _level = other._level;
  _name = other._name;
  _type = other._type;
}

bool tensor_index::operator == (const tensor_index& A) const{
  return (_size==A._size)&&(_name==A._name)&&(_type==A._type)&&(_level==A._level);
}

bool tensor_index::operator == (tensor_index& A){
  return (_size==A._size)&&(_name==A._name)&&(_type==A._type)&&(_level==A._level);
}

inline bool tensor_index::similar(const tensor_index& A){
  return (_size==A._size)&&(_type==A._type)&&(_level!=A._level);
}

inline bool tensor_index::similar(tensor_index& A){
  return (_size==A._size)&&(_type==A._type)&&(_level!=A._level);
}

inline bool tensor_index::similar(const tensor_index& A, index_type type){
  return (_size==A._size)&&(_type==A._type)&&(_type==type)&&(_level!=A._level);
}

inline bool tensor_index::similar(tensor_index& A, index_type type){
  return (_size==A._size)&&(_type==A._type)&&(_type==type)&&(_level!=A._level);
}

#endif
