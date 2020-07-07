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


#ifndef INDEX_CLASS_FOR_DENSE_TENSOR_HEADER
#define INDEX_CLASS_FOR_DENSE_TENSOR_HEADER

#include "../util/types_and_headers.h"

class dtensor_index{
public:
  dtensor_index();
  dtensor_index(int64_t size);
  dtensor_index(int64_t size, string name);
  dtensor_index(int64_t size, string name, index_type type);
  dtensor_index(int64_t size, string name, index_type type, unsigned level);
  dtensor_index(const dtensor_index& other);

  ~dtensor_index(){}

  inline void setTag(){_tag = to_string(_size)+" "+_name+" "+to_string(_type)+" ";}
  inline void rename(string name){_name=name; setTag();}
  inline void resize(int64_t size){_size=size; setTag();}

  dtensor_index operator = (const dtensor_index& other);
  bool operator == (const dtensor_index& A) const;
  bool operator != (const dtensor_index& A) const;
  bool operator <  (const dtensor_index& A) const;

  inline int64_t size(){return _size;}
  inline unsigned level(){return _level;}
  inline string name(){return _name;}
  inline index_type type(){return _type;}
  inline string tag(){return _tag+to_string(_level);}

  inline bool similar(const dtensor_index& A){
    return (_size==A._size)&&(_name==A._name)&&(_type==A._type)&&((_level==A._level+1)||(_level+1==A._level));
  }

  inline bool similar(dtensor_index& A){
    return (_size==A._size)&&(_name==A._name)&&(_type==A._type)&&((_level==A._level+1)||(_level+1==A._level));
  }

  inline bool similar(const dtensor_index& A, index_type type){
    return (_size==A._size)&&(_name==A._name)&&(_type==A._type || type==All)&&(_type==type || type==All)&&((_level==A._level+1)||(_level+1==A._level));
  }

  inline bool similar(dtensor_index& A, index_type type){
    return (_size==A._size)&&(_name==A._name)&&(_type==A._type || type==All)&&(_type==type || type==All)&&((_level==A._level+1)||(_level+1==A._level));
  }

  inline void prime(int inc=1){_level+=inc;}
  inline void primeLink(int inc=1){if(_type==Link) _level+=inc;}
  inline void primeSite(int inc=1){if(_type==Site) _level+=inc;}
  inline void mapPrime(unsigned from, unsigned to){if(_level==from) _level=to;}
  inline void mapPrime(unsigned from, unsigned to, index_type type){if(_level==from&&(_type==type||type==All)) _level=to;}
  inline void noPrime(index_type type=All){if(type==All || _type==type) _level=0;}

private:
  int64_t _size;
  unsigned _level;
  index_type _type;
  string _name;
  string _tag; // to make it hashable
};

dtensor_index noPrime(dtensor_index& index);

#endif
