#ifndef QN_TYPE_HEADER
#define QN_TYPE_HEADER

// describes a system with 1 abelian quantum number
class QN_t_1 {

  using qtype = int;
  public:
  QN_t_1(): qn(0) {};
  QN_t_1(qtype qn_in): qn(qn_in) {};
  //QN_t_1(std::initialize_list<qtype> qn_list){ qn = qn_list; }

  inline QN_t_1 operator *  (const QN_t_1& A) const{ return {qn*A.qn}; }
  inline QN_t_1 operator +  (const QN_t_1& A) const{ return {qn+A.qn}; }
  //inline QN_t_1 operator +  (const QN_t_1 A) const{  return {qn+A.qn}; }
  inline QN_t_1 operator -  (const QN_t_1& A) const{ return {qn-A.qn};};
  inline QN_t_1 operator - (){ return {-qn};};
  QN_t_1& operator += (const QN_t_1& A){ qn+=A.qn; return *this; }
  QN_t_1& operator -= (const QN_t_1& A){ qn-=A.qn; return *this; }
  inline bool   operator == (const QN_t_1& rhs)const{ return qn == rhs.qn;    }
  inline bool   operator != (const QN_t_1& rhs)const{ return !(*this == rhs); }
  inline bool   operator <  (const QN_t_1& rhs)const{ return qn < rhs.qn;     }

  /*qtype operator ()() const{
   *     return qn;
   *       } */

  friend std::ostream& operator<<(std::ostream &os, const QN_t_1 &A){
    os << "(" << A.qn << ")";
    return os;
  }

  //Note: illegal to inject to_string customization into std namespace because it isn't
  friend std::string to_string(const QN_t_1& A){ std::string s = std::to_string(A.qn); return s;  }

  qtype qn;
};
namespace std {
  template <> struct hash<QN_t_1>{
    size_t operator()(const QN_t_1& A) const noexcept{
      return hash<int>{}(A.qn);
    }
  };
}
#endif
