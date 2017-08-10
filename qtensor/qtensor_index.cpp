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

inline bool qtensor_index::similar(const qtensor_index& A){
  if(_arrow != A._arrow &&
    _tag    == A._tag   &&
    ((_level== A._level+1)||(_level+1== A._level))
  ){
    return true;
  }else{
    return false;
  }
}

inline bool qtensor_index::similar(qtensor_index& A){
  if(_arrow != A._arrow &&
    _tag    == A._tag   &&
    ((_level== A._level+1)||(_level+1== A._level))
  ){
    return true;
  }else{
    return false;
  }
}

bool qtensor_index::operator == (const qtensor_index& A) const{
  return (_tag==A._tag && _arrow==A._arrow && _level==A._level);
}

bool qtensor_index::operator != (const qtensor_index& A) const{
  return (_tag!=A._tag || _arrow!=A._arrow || _level!=A._level);
}

bool qtensor_index::operator <  (const qtensor_index& A) const{
  return (_tag<A._tag && _arrow<A._arrow && _level<A._level);
}

#endif
