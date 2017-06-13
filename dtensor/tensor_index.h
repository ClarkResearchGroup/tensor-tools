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
  bool operator == (tensor_index& A);
  inline unsigned size(){return _size;}
  inline unsigned level(){return _level;}
  inline string name(){return _name;}
  inline index_type type(){return _type;}

  inline bool similar(const tensor_index& A);
  inline bool similar(tensor_index& A);
  inline bool similar(const tensor_index& A, index_type type);
  inline bool similar(tensor_index& A, index_type type);

  inline void prime(int inc=1){_level+=inc;}
  inline void primeLink(){if(_type==Link) ++_level;}
  inline void primeSite(){if(_type==Site) ++_level;}
  inline void mapPrime(unsigned from, unsigned to){if(_level==from) _level=to;}
  inline void mapPrime(unsigned from, unsigned to, index_type type){if(_level==from&&_type==type) _level=to;}

  unsigned _size;
  unsigned _level;
  index_type _type;
  string _name;
};

#endif
