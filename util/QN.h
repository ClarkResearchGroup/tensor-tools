#ifndef QN_TYPE_HEADER
#define QN_TYPE_HEADER
#include <valarray>

// describes a system with 1 abelian quantum number
class QN_t_1 {

  //using qtype = int;
  using qtype = std::valarray<int>;
  public:
  QN_t_1(): qn(0) {};
  QN_t_1(qtype qn_in): qn(qn_in) {};
  QN_t_1(int qn_in):qn({qn_in}) {};
  QN_t_1(std::initializer_list<int> qn_list): qn(qn_list){};

  inline QN_t_1 operator *  (const QN_t_1& A) const{ return {qn*A.qn}; }
  inline QN_t_1 operator +  (const QN_t_1& A) const{ return {qn+A.qn}; }
  //inline QN_t_1 operator +  (const QN_t_1 A) const{  return {qn+A.qn}; }
  inline QN_t_1 operator -  (const QN_t_1& A) const{ return {qn-A.qn};};
  inline QN_t_1 operator - (){ return {-qn};};
  QN_t_1& operator += (const QN_t_1& A){ qn+=A.qn; return *this; }
  QN_t_1& operator -= (const QN_t_1& A){ qn-=A.qn; return *this; }
  inline bool operator == (const QN_t_1& rhs)const{ 
    assert(rhs.qn.size()==qn.size());
    return (qn == rhs.qn).sum()==qn.size();    
  }
  inline bool   operator != (const QN_t_1& rhs)const{ return !(*this == rhs); }
  inline bool   operator <  (const QN_t_1& rhs)const{ 
    //return qn < rhs.qn;     
    return qn[0] < rhs.qn[0];  
  }

  /*qtype operator ()() const{
   *     return qn;
   *       } */

  friend std::ostream& operator<<(std::ostream &os, const QN_t_1 &A){
    os << "(";
   for(auto c: A.qn) 
     os<< c <<",";
    os<< ")";
    return os;
  }

  //Note: illegal to inject to_string customization into std namespace because it isn't
  friend std::string to_string(const QN_t_1& A){ 
    std::string s;
    for(auto c: A.qn) s+=std::to_string(c)+"|"; 
    return s;  
  }

  qtype qn;
};
namespace std {
  template <> struct hash<QN_t_1>{
    size_t operator()(const QN_t_1& A) const noexcept{
      return hash<std::string>{}(to_string(A));
    }
  };
}
// describes a system with 2 abelian quantum number
class QN_t_2 {

  //using qtype = int;
  using qtype = std::valarray<int>;
  public:
  QN_t_2(): qn({0,0}) {};
  QN_t_2(int qn_in):qn({qn_in,0}) {assert(qn_in==0);};
  QN_t_2(qtype qn_in): qn(qn_in) {};
  QN_t_2(std::initializer_list<int> qn_list): qn(qn_list){};

  inline QN_t_2 operator *  (const QN_t_2& A) const{ return {qn*A.qn}; }
  inline QN_t_2 operator +  (const QN_t_2& A) const{ return {qn+A.qn}; }
  inline QN_t_2 operator -  (const QN_t_2& A) const{ return {qn-A.qn};};
  inline QN_t_2 operator - (){ return {-qn};};
  QN_t_2& operator += (const QN_t_2& A){ qn+=A.qn; return *this; }
  QN_t_2& operator -= (const QN_t_2& A){ qn-=A.qn; return *this; }
  inline bool operator == (const QN_t_2& rhs)const{ 
    assert(rhs.qn.size()==qn.size());
    return (qn[0] == rhs.qn[0]) && (qn[1]==rhs.qn[1]);
  }
  inline bool   operator != (const QN_t_2& rhs)const{ return !(*this == rhs); }
  inline bool   operator <  (const QN_t_2& rhs)const{ 
    //return qn < rhs.qn;     
    //return qn[0] < rhs.qn[0];  
    return std::make_pair(qn[0],qn[1]) < std::make_pair(rhs.qn[0],rhs.qn[1]);
  }

  /*qtype operator ()() const{
   *     return qn;
   *       } */

  friend std::ostream& operator<<(std::ostream &os, const QN_t_2 &A){
    os << "(";
   for(auto c: A.qn) 
     os<< c <<",";
    os<< ")";
    return os;
  }

  //Note: illegal to inject to_string customization into std namespace because it isn't
  friend std::string to_string(const QN_t_2 A){ 
    std::string s;
    for(auto c: A.qn) s+=std::to_string(c)+"|"; 
    return s;  
  }

  qtype qn;
};
namespace std {
  template <> struct hash<QN_t_2>{
    size_t operator()(const QN_t_2& A) const noexcept{
      return hash<std::string>{}(to_string(A));
    }
  };
}
#endif
