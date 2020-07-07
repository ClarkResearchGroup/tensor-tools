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

class tensor_index{
public:
  tensor_index(){}
  tensor_index(unsigned size);
  tensor_index(unsigned size, string name);
  tensor_index(unsigned size, string name, index_type type);
  tensor_index(unsigned size, string name, index_type type, unsigned level);
  tensor_index(const tensor_index& other);

  ~tensor_index(){}

  bool operator == (const tensor_index& A) const;
  bool operator != (const tensor_index& A) const;
  bool operator <  (const tensor_index& A) const;

  inline unsigned size(){return _size;}
  inline unsigned level(){return _level;}
  inline string name(){return _name;}
  inline index_type type(){return _type;}

  inline bool similar(const tensor_index& A);
  inline bool similar(tensor_index& A);
  inline bool similar(const tensor_index& A, index_type type);
  inline bool similar(tensor_index& A, index_type type);

  inline void prime(int inc=1){_level+=inc; _tag = to_string(_size)+" "+_name+" "+to_string(_type)+" "+to_string(_level);}
  inline void primeLink(){if(_type==Link) ++_level; _tag = to_string(_size)+" "+_name+" "+to_string(_type)+" "+to_string(_level);}
  inline void primeSite(){if(_type==Site) ++_level; _tag = to_string(_size)+" "+_name+" "+to_string(_type)+" "+to_string(_level);}
  inline void mapPrime(unsigned from, unsigned to){if(_level==from) _level=to; _tag = to_string(_size)+" "+_name+" "+to_string(_type)+" "+to_string(_level);}
  inline void mapPrime(unsigned from, unsigned to, index_type type){if(_level==from&&_type==type) _level=to; _tag = to_string(_size)+" "+_name+" "+to_string(_type)+" "+to_string(_level);}

  unsigned _size;
  unsigned _level;
  index_type _type;
  string _name;

  string _tag; // to make it hashable
};

#endif
