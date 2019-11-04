#ifndef ALL_MPS_RELATED_HEADERS
#define ALL_MPS_RELATED_HEADERS

#include "tt.h"
#include "qtt.h"
#include "qstt.h"
#include "observables.h"
#include "mps_mpo_methods.h"
namespace ezh5{
  template<> Node& Node::operator = (QN_t_1 val){
    (*this) = val.qn[0];
    return *this;
  }
  template<> Node& Node::operator >> (QN_t_1& val){
    int temp;
    (*this) >> temp;
    val.qn = {temp};
    return *this;
  }
  template<> Node& Node::operator = (std::vector<QN_t_1>& vec){
    vector<int> temp;temp.reserve(vec.size()); for(auto& c: vec) temp.push_back(c.qn[0]);
    (*this) = temp;
    return *this;
  };
  template<> Node& Node::operator >> (std::vector<QN_t_1>& vec){
    vector<int> temp; (*this) >> temp;
    for(auto c: temp) vec.push_back({c});
    return *this;
  };
  //----------------------------------------------------------------
  template<> Node& Node::operator = (QN_t_2 val){
    vector<int> temp;temp.reserve(val.qn.size()); for(auto& c: val.qn) temp.push_back(c);
    (*this) = temp;
    return *this;
  }
  template<> Node& Node::operator >> (QN_t_2& val){
    vector<int> temp; (*this) >> temp;
    val.qn = valarray<int>(temp.data(),temp.size());
    return *this;
  }
  template<> Node& Node::operator = (std::vector<QN_t_2>& vec){
    vector<int> temp;temp.reserve(vec.size()); 
    for(auto& v: vec) 
      for(auto c: v.qn)
        temp.push_back(c);
    (*this) = temp;
    return *this;
  };
  template<> Node& Node::operator >> (std::vector<QN_t_2>& vec){
    vector<int> temp; (*this) >> temp;
    assert(temp.size()%2 == 0);
    for(size_t i = 0; i<temp.size(); i+=2) 
      vec.push_back({temp[i],temp[i+1]});
    return *this;
  };
}

#endif
