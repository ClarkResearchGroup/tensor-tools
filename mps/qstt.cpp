#ifndef QUANTUM_NUMBERED_SPARSE_TENSORTRAIN_CLASS
#define QUANTUM_NUMBERED_SPARSE_TENSORTRAIN_CLASS

#include "qstt.h"

//---------------------------------------------------------------------------
// Constructors
template <typename T, unsigned N>
qsTensorTrain<T, N>::qsTensorTrain() {
  tensors_allocated = false;
  length  = 0;
  phy_dim = 0;
  center  = -1;
  totalQ  = 0;
  _id = unsigned(1e9*thread_safe_random_double()+1e7*thread_safe_random_double()+1e5*thread_safe_random_double()+1e3*thread_safe_random_double());
}
template qsTensorTrain<double, 1>::qsTensorTrain();
template qsTensorTrain<double, 2>::qsTensorTrain();
template qsTensorTrain<std::complex<double>, 1>::qsTensorTrain();
template qsTensorTrain<std::complex<double>, 2>::qsTensorTrain();


template <typename T, unsigned N>
qsTensorTrain<T, N>::qsTensorTrain(abstract_sites* s, int Q){
  _id = unsigned(1e9*thread_safe_random_double()+1e7*thread_safe_random_double()+1e5*thread_safe_random_double()+1e3*thread_safe_random_double());
  tensors_allocated = false;
  setLength((*s).N());
  setPhysicalDim((*s).phy_dim());
  for (size_t i = 0; i < length+1; i++) {
    bond_dims[i] = 1;
  }
  phy_qn = (*s).phy_qn();
  totalQ = Q;
  allocateTensors();
  center = -1;
}
template qsTensorTrain<double, 1>::qsTensorTrain(abstract_sites* s, int Q);
template qsTensorTrain<double, 2>::qsTensorTrain(abstract_sites* s, int Q);
template qsTensorTrain<std::complex<double>, 1>::qsTensorTrain(abstract_sites* s, int Q);
template qsTensorTrain<std::complex<double>, 2>::qsTensorTrain(abstract_sites* s, int Q);


template <typename T, unsigned N>
qsTensorTrain<T, N>::qsTensorTrain(abstract_sites* s, str_vec product_string){
  _id = unsigned(1e9*thread_safe_random_double()+1e7*thread_safe_random_double()+1e5*thread_safe_random_double()+1e3*thread_safe_random_double());
  tensors_allocated = false;
  setLength((*s).N());
  setPhysicalDim((*s).phy_dim());
  for (size_t i = 0; i < length+1; i++) {
    bond_dims[i] = 1;
  }
  phy_qn = (*s).phy_qn();
  uint_vec product_state = (*s).product_state(product_string);
  totalQ = 0;
  assert(product_state.size()==length);
  for (size_t i = 0; i < length; i++) {
    totalQ += phy_qn[product_state[i]];
  }
  allocateTensors(product_state.data());
  center = -1;
}
template qsTensorTrain<double, 1>::qsTensorTrain(abstract_sites* s, str_vec product_string);
template qsTensorTrain<double, 2>::qsTensorTrain(abstract_sites* s, str_vec product_string);
template qsTensorTrain<std::complex<double>, 1>::qsTensorTrain(abstract_sites* s, str_vec product_string);
template qsTensorTrain<std::complex<double>, 2>::qsTensorTrain(abstract_sites* s, str_vec product_string);


template <typename T, unsigned N>
qsTensorTrain<T, N>::qsTensorTrain(unsigned L, unsigned pD, vector<int>& phy_qn, int Q){
  _id = unsigned(1e9*thread_safe_random_double()+1e7*thread_safe_random_double()+1e5*thread_safe_random_double()+1e3*thread_safe_random_double());
  tensors_allocated = false;
  setLength(L);
  setPhysicalDim(pD);
  for (size_t i = 0; i < length+1; i++) {
    bond_dims[i] = 1;
  }
  phy_qn = phy_qn;
  totalQ = Q;
  allocateTensors();
  center = -1;
}
template qsTensorTrain<double, 1>::qsTensorTrain(unsigned L, unsigned pD, vector<int>& phy_qn, int totalQ);
template qsTensorTrain<double, 2>::qsTensorTrain(unsigned L, unsigned pD, vector<int>& phy_qn, int totalQ);
template qsTensorTrain<std::complex<double>, 1>::qsTensorTrain(unsigned L, unsigned pD, vector<int>& phy_qn, int totalQ);
template qsTensorTrain<std::complex<double>, 2>::qsTensorTrain(unsigned L, unsigned pD, vector<int>& phy_qn, int totalQ);


template <typename T, unsigned N>
qsTensorTrain<T, N>::qsTensorTrain(unsigned L, unsigned pD, vector<int>& phy_qn, uint_vec product_state){
  _id = unsigned(1e9*thread_safe_random_double()+1e7*thread_safe_random_double()+1e5*thread_safe_random_double()+1e3*thread_safe_random_double());
  tensors_allocated = false;
  setLength(L);
  setPhysicalDim(pD);
  for (size_t i = 0; i < length+1; i++) {
    bond_dims[i] = 1;
  }
  phy_qn = phy_qn;
  totalQ = 0;
  assert(product_state.size()==length);
  for (size_t i = 0; i < length; i++) {
    totalQ += phy_qn[product_state[i]];
  }
  allocateTensors(product_state.data());
  center = -1;
}
template qsTensorTrain<double, 1>::qsTensorTrain(unsigned L, unsigned pD, vector<int>& phy_qn, uint_vec product_state);
template qsTensorTrain<double, 2>::qsTensorTrain(unsigned L, unsigned pD, vector<int>& phy_qn, uint_vec product_state);
template qsTensorTrain<std::complex<double>, 1>::qsTensorTrain(unsigned L, unsigned pD, vector<int>& phy_qn, uint_vec product_state);
template qsTensorTrain<std::complex<double>, 2>::qsTensorTrain(unsigned L, unsigned pD, vector<int>& phy_qn, uint_vec product_state);


template <typename T, unsigned N>
qsTensorTrain<T, N>::qsTensorTrain(const qsTensorTrain<T, N>& other){
  tensors_allocated = false;
  if(other.tensors_allocated){
    _id = other._id;
    setLength(other.length);
    setPhysicalDim(other.phy_dim);
    bond_dims = other.bond_dims;
    totalQ = other.totalQ;
    phy_qn = other.phy_qn;
    center = other.center;
    A.clear();
    A = other.A;
    tensors_allocated = true;
  }
}
template qsTensorTrain<double, 1>::qsTensorTrain(const qsTensorTrain<double, 1>& other);
template qsTensorTrain<double, 2>::qsTensorTrain(const qsTensorTrain<double, 2>& other);
template qsTensorTrain<std::complex<double>, 1>::qsTensorTrain(const qsTensorTrain< std::complex<double>, 1>& other);
template qsTensorTrain<std::complex<double>, 2>::qsTensorTrain(const qsTensorTrain< std::complex<double>, 2>& other);


template <typename T, unsigned N>
qsTensorTrain<T, N>::qsTensorTrain(qsTensorTrain<T, N>&& other){
  tensors_allocated = false;
  if(other.tensors_allocated){
    _id = other._id;
    setLength(other.length);
    setPhysicalDim(other.phy_dim);
    bond_dims = std::move(other.bond_dims);
    totalQ = other.totalQ;
    phy_qn = std::move(other.phy_qn);
    center = other.center;
    A.clear();
    A = std::move(other.A);
    tensors_allocated = true;
    other.tensors_allocated = false;
  }
}
template qsTensorTrain<double, 1>::qsTensorTrain(qsTensorTrain<double, 1>&& other);
template qsTensorTrain<double, 2>::qsTensorTrain(qsTensorTrain<double, 2>&& other);
template qsTensorTrain<std::complex<double>, 1>::qsTensorTrain(qsTensorTrain< std::complex<double>, 1>&& other);
template qsTensorTrain<std::complex<double>, 2>::qsTensorTrain(qsTensorTrain< std::complex<double>, 2>&& other);

template <typename T, unsigned N>
qsTensorTrain<T, N>::qsTensorTrain(const qTensorTrain<T, N>& other){
  tensors_allocated = false;
  if(other.tensors_allocated){
    _id = other._id;
    setLength(other.length);
    setPhysicalDim(other.phy_dim);
    bond_dims = other.bond_dims;
    totalQ = other.totalQ;
    phy_qn = other.phy_qn;
    center = other.center;
    A.clear();
    //A = other.A;
    A.resize(other.A.size());
    for(int i=0;i<A.size();i++)
      A[i]= other.A[i];
    tensors_allocated = true;
  }
}
template qsTensorTrain<double, 1>::qsTensorTrain(const qTensorTrain<double, 1>& other);
template qsTensorTrain<double, 2>::qsTensorTrain(const qTensorTrain<double, 2>& other);
template qsTensorTrain<std::complex<double>, 1>::qsTensorTrain(const qTensorTrain< std::complex<double>, 1>& other);
template qsTensorTrain<std::complex<double>, 2>::qsTensorTrain(const qTensorTrain< std::complex<double>, 2>& other);


template <typename T, unsigned N>
qsTensorTrain<T, N>::qsTensorTrain(qTensorTrain<T, N>&& other){
  tensors_allocated = false;
  if(other.tensors_allocated){
    _id = other._id;
    setLength(other.length);
    setPhysicalDim(other.phy_dim);
    bond_dims = std::move(other.bond_dims);
    totalQ = other.totalQ;
    phy_qn = std::move(other.phy_qn);
    center = other.center;
    A.clear();
    //A = std::move(other.A);
    A.resize(other.A.size());
    for(int i=0;i<A.size();i++)
      A[i] = other.A[i];
    tensors_allocated = true;
    other.tensors_allocated = false;
  }
}
template qsTensorTrain<double, 1>::qsTensorTrain(qTensorTrain<double, 1>&& other);
template qsTensorTrain<double, 2>::qsTensorTrain(qTensorTrain<double, 2>&& other);
template qsTensorTrain<std::complex<double>, 1>::qsTensorTrain(qTensorTrain< std::complex<double>, 1>&& other);
template qsTensorTrain<std::complex<double>, 2>::qsTensorTrain(qTensorTrain< std::complex<double>, 2>&& other);
//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
// Set shapes
template <typename T, unsigned N>
void qsTensorTrain<T, N>::setLength(int L) {
  assert(L>2);
  assert(!tensors_allocated);
  length = L;
  bond_dims.resize(length+1);
}
template void qsTensorTrain<double, 1>::setLength(int L);
template void qsTensorTrain<double, 2>::setLength(int L);
template void qsTensorTrain<std::complex<double>, 1>::setLength(int L);
template void qsTensorTrain<std::complex<double>, 2>::setLength(int L);


template <typename T, unsigned N>
void qsTensorTrain<T, N>::setPhysicalDim(int s){
  assert(s>0);
  assert(!tensors_allocated);
  phy_dim = s;
}
template void qsTensorTrain<double, 1>::setPhysicalDim(int s);
template void qsTensorTrain<double, 2>::setPhysicalDim(int s);
template void qsTensorTrain<std::complex<double>, 1>::setPhysicalDim(int s);
template void qsTensorTrain<std::complex<double>, 2>::setPhysicalDim(int s);
//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
// Tensor data management
template <typename T, unsigned N>
void qsTensorTrain<T, N>::allocateTensors(unsigned* product_state){
  if(!tensors_allocated) {
    string Link_name_pref = "ID"+to_string(_id)+"Link";
    string Site_name_pref = "Site";
    if(N==1){
      // MPS
      // Get legal Link QN
      vector< set<int> > Link_QN(length+1);
      Link_QN[0].insert(0);
      Link_QN[length].insert(totalQ);
      for (size_t j = 1; j < length; j++) {
        for (set<int>::iterator it=Link_QN[j-1].begin(); it!=Link_QN[j-1].end(); ++it) {
          for (size_t k = 0; k < phy_qn.size(); k++) {
            Link_QN[j].insert(*it + phy_qn[k]);
          }
        }
      }
      for (size_t j = length-1; j > 0; j--) {
        vector<int> to_remove;
        for (set<int>::iterator it=Link_QN[j].begin(); it!=Link_QN[j].end(); ++it) {
          bool valid = false;
          for (size_t k = 0; k < phy_qn.size(); k++) {
            if(Link_QN[j+1].find(*it + phy_qn[k]) != Link_QN[j+1].end()){
              valid = true;
              break;
            }
          }
          if(!valid){
            to_remove.push_back(*it);
          }
        }
        for (size_t k = 0; k < to_remove.size(); k++) {
          Link_QN[j].erase(to_remove[k]);
        }
      }
      // Init tensors
      for (size_t i = 0; i < length; i++) {
        string left_link_name  = Link_name_pref+to_string(i);
        string right_link_name = Link_name_pref+to_string(i+1);
        string site_name       = Site_name_pref+to_string(i);
        A.push_back(
          std::move(
            qstensor<T>(
              {Inward, Inward, Outward},
              {left_link_name, site_name, right_link_name},
              {Link, Site, Link},
              {0,0,0})
          )
        );
        // Set QNs
        vector<int> left_link_qn(Link_QN[i].begin(), Link_QN[i].end());
        vector<int> right_link_qn(Link_QN[i+1].begin(), Link_QN[i+1].end());
        // Left Link
        for (size_t j = 0; j < left_link_qn.size(); j++) {
          A[i].addQNtoIndex(0, std::make_pair(left_link_qn[j], 1));
        }
        // Site
        for (size_t j = 0; j < phy_qn.size(); j++) {
          A[i].addQNtoIndex(1, std::make_pair(phy_qn[j], 1));
        }
        // Right Link
        for (size_t j = 0; j < right_link_qn.size(); j++) {
          A[i].addQNtoIndex(2, std::make_pair(right_link_qn[j], 1));
        }
        A[i].initBlock();
        //A[i].setRandom();
      }
      // Set product state
      if(product_state!=nullptr){
        for (size_t i = 0; i < length; i++) {
          string ind = A[i].getIndices();
          for (size_t j = 0; j < A[i]._block.size(); j++) {
            if(A[i].block_index_qi[j][1] == product_state[i]){
              A[i]._block[j][ind.c_str()] = 1.;
            }else{
              //A[i]._block[j][ind.c_str()] = 0;//std::numeric_limits<double>::epsilon(); //TODO: check for complex
            }
          }
        }
      }
    }
    if(N==2){
      // MPO is always init to identity
      for (size_t i = 0; i < length; i++) {
        string left_link_name  = Link_name_pref+to_string(i);
        string right_link_name = Link_name_pref+to_string(i+1);
        string site_name       = Site_name_pref+to_string(i);
        A.push_back(
          std::move(
            qstensor<T>(
              {Inward, Outward, Inward, Outward},
              {left_link_name, site_name, site_name, right_link_name},
              {Link, Site, Site, Link},
              {0,0,1,0})
          )
        );
        // Set QNs
        // Left Link
        A[i].addQNtoIndex(0, std::make_pair(0, 1));
        // Site
        for (size_t j = 0; j < phy_qn.size(); j++) {
          A[i].addQNtoIndex(1, std::make_pair(phy_qn[j], 1));
        }
        // Site
        for (size_t j = 0; j < phy_qn.size(); j++) {
          A[i].addQNtoIndex(2, std::make_pair(phy_qn[j], 1));
        }
        // Right Link
        A[i].addQNtoIndex(3, std::make_pair(0, 1));
        // Set up blocks
        A[i].initBlock();
        //A[i].setZero();
        // make all diagonal blocks equal
        /*for (size_t j = 0; j < A[i].block.size(); j++) {
          for (size_t k = 0; k < A[i].block[j].size(); k++) {
            A[i].block[j][k] = 1.0;
          }
        }*/
        A[i].setOne();
      }
    }
    tensors_allocated = true;
  }
}
// template <typename T, unsigned N>
// void qsTensorTrain<T, N>::allocateTensors(unsigned* product_state){
//   if(!tensors_allocated) {
//     string Link_name_pref = "ID"+to_string(_id)+"Link";
//     string Site_name_pref = "Site";
//     if(N==1){
//       // MPS
//       if(product_state==nullptr){
//         // Get legal Link QN
//         vector< set<int> > Link_QN(length+1);
//         Link_QN[0].insert(0);
//         Link_QN[length].insert(totalQ);
//         for (size_t j = 1; j < length; j++) {
//           for (set<int>::iterator it=Link_QN[j-1].begin(); it!=Link_QN[j-1].end(); ++it) {
//             for (size_t k = 0; k < phy_qn.size(); k++) {
//               Link_QN[j].insert(*it + phy_qn[k]);
//             }
//           }
//         }
//         for (size_t j = length-1; j > 0; j--) {
//           vector<int> to_remove;
//           for (set<int>::iterator it=Link_QN[j].begin(); it!=Link_QN[j].end(); ++it) {
//             bool valid = false;
//             for (size_t k = 0; k < phy_qn.size(); k++) {
//               if(Link_QN[j+1].find(*it + phy_qn[k]) != Link_QN[j+1].end()){
//                 valid = true;
//                 break;
//               }
//             }
//             if(!valid){
//               to_remove.push_back(*it);
//             }
//           }
//           for (size_t k = 0; k < to_remove.size(); k++) {
//             Link_QN[j].erase(to_remove[k]);
//           }
//         }
//         // Init tensors
//         for (size_t i = 0; i < length; i++) {
//           string left_link_name  = Link_name_pref+to_string(i);
//           string right_link_name = Link_name_pref+to_string(i+1);
//           string site_name       = Site_name_pref+to_string(i);
//           A.push_back(
//             std::move(
//               qstensor<T>(
//                 {Inward, Inward, Outward},
//                 {left_link_name, site_name, right_link_name},
//                 {Link, Site, Link},
//                 {0,0,0})
//             )
//           );
//           // Set QNs
//           vector<int> left_link_qn(Link_QN[i].begin(), Link_QN[i].end());
//           vector<int> right_link_qn(Link_QN[i+1].begin(), Link_QN[i+1].end());
//           // Left Link
//           for (size_t j = 0; j < left_link_qn.size(); j++) {
//             A[i].addQNtoIndex(0, std::make_pair(left_link_qn[j], 1));
//           }
//           // Site
//           for (size_t j = 0; j < phy_qn.size(); j++) {
//             A[i].addQNtoIndex(1, std::make_pair(phy_qn[j], 1));
//           }
//           // Right Link
//           for (size_t j = 0; j < right_link_qn.size(); j++) {
//             A[i].addQNtoIndex(2, std::make_pair(right_link_qn[j], 1));
//           }
//           A[i].initBlock();
//           A[i].setRandom();
//         }
//       }else{
//         int_vec Link_QN(length+1);
//         Link_QN[0] = 0;
//         for (size_t j = 0; j < length; j++) {
//           for (size_t k = 0; k < phy_qn.size(); k++) {
//             if(k==product_state[j]) Link_QN[j+1] = Link_QN[j] + phy_qn[k];
//           }
//         }
//         // Init tensors
//         for (size_t i = 0; i < length; i++) {
//           string left_link_name  = Link_name_pref+to_string(i);
//           string right_link_name = Link_name_pref+to_string(i+1);
//           string site_name       = Site_name_pref+to_string(i);
//           A.push_back(
//             std::move(
//               qstensor<T>(
//                 {Inward, Inward, Outward},
//                 {left_link_name, site_name, right_link_name},
//                 {Link, Site, Link},
//                 {0,0,0})
//             )
//           );
//           // Left Link
//           A[i].addQNtoIndex(0, std::make_pair(Link_QN[i], 1));
//           // Site
//           for (size_t j = 0; j < phy_qn.size(); j++) {
//             A[i].addQNtoIndex(1, std::make_pair(phy_qn[j], 1));
//           }
//           // Right Link
//           A[i].addQNtoIndex(2, std::make_pair(Link_QN[i+1], 1));
//           A[i].initBlock();
//           A[i].setOne();
//         }
//       }
//     }
//     if(N==2){
//       // MPO is always init to identity
//       for (size_t i = 0; i < length; i++) {
//         string left_link_name  = Link_name_pref+to_string(i);
//         string right_link_name = Link_name_pref+to_string(i+1);
//         string site_name       = Site_name_pref+to_string(i);
//         A.push_back(
//           std::move(
//             qstensor<T>(
//               {Inward, Outward, Inward, Outward},
//               {left_link_name, site_name, site_name, right_link_name},
//               {Link, Site, Site, Link},
//               {0,0,1,0})
//           )
//         );
//         // Set QNs
//         // Left Link
//         A[i].addQNtoIndex(0, std::make_pair(0, 1));
//         // Site
//         for (size_t j = 0; j < phy_qn.size(); j++) {
//           A[i].addQNtoIndex(1, std::make_pair(phy_qn[j], 1));
//         }
//         // Site
//         for (size_t j = 0; j < phy_qn.size(); j++) {
//           A[i].addQNtoIndex(2, std::make_pair(phy_qn[j], 1));
//         }
//         // Right Link
//         A[i].addQNtoIndex(3, std::make_pair(0, 1));
//         // Set up blocks
//         A[i].initBlock();
//         A[i].setZero();
//         // make all diagonal blocks equal
//         for (size_t j = 0; j < A[i].block.size(); j++) {
//           for (size_t k = 0; k < A[i].block[j].size(); k++) {
//             A[i].block[j].data()[k] = 1.0;
//           }
//         }
//       }
//     }
//     tensors_allocated = true;
//   }
// }
template void qsTensorTrain<double, 1>::allocateTensors(unsigned* product_state);
template void qsTensorTrain<double, 2>::allocateTensors(unsigned* product_state);
template void qsTensorTrain<std::complex<double>, 1>::allocateTensors(unsigned* product_state);
template void qsTensorTrain<std::complex<double>, 2>::allocateTensors(unsigned* product_state);


template <typename T, unsigned N>
void qsTensorTrain<T, N>::freeTensors(){
  if(tensors_allocated) {
    A.clear();
    tensors_allocated = false;
  }
}
template void qsTensorTrain<double, 1>::freeTensors();
template void qsTensorTrain<double, 2>::freeTensors();
template void qsTensorTrain<std::complex<double>, 1>::freeTensors();
template void qsTensorTrain<std::complex<double>, 2>::freeTensors();


//---------------------------------------------------------------------------
// Set tensors values
template <typename T, unsigned N>
void qsTensorTrain<T, N>::setZero(){
  if(tensors_allocated) {
    for (size_t i = 0; i < length; i++) {
      A[i].setZero();
    }
  }
}
template void qsTensorTrain<double, 1>::setZero();
template void qsTensorTrain<double, 2>::setZero();
template void qsTensorTrain<std::complex<double>, 1>::setZero();
template void qsTensorTrain<std::complex<double>, 2>::setZero();


template <typename T, unsigned N>
void qsTensorTrain<T, N>::setRandom(){
  if(tensors_allocated) {
    for (size_t i = 0; i < length; i++) {
      A[i].setRandom();
    }
  }
}
template void qsTensorTrain<double, 1>::setRandom();
template void qsTensorTrain<double, 2>::setRandom();
template void qsTensorTrain<std::complex<double>, 1>::setRandom();
template void qsTensorTrain<std::complex<double>, 2>::setRandom();


//---------------------------------------------------------------------------
// Print qsTensorTrain information
template <typename T, unsigned N>
void qsTensorTrain<T, N>::print(int level){
  pout<<"---------------------------------------------------"<<'\n';
  pout << "length = " << length << ",  phy_dim = " << phy_dim << ", totalQ = " << totalQ << ", center = " << center << ", id = " << _id << '\n';
  pout << "bond_dims vector = " << " ";
  for(auto v : bond_dims) pout << v << " ";
  pout<<std::endl;
  if(level>=1){
    pout << "Information of individual tensors:" << '\n';
    if(tensors_allocated){
      for (size_t i = 0; i < length; i++) {
        A[i].print(std::max(level-1,0));
      }
    }
  }
  pout<<"---------------------------------------------------"<<'\n';
}
template void qsTensorTrain<double, 1>::print(int level);
template void qsTensorTrain<double, 2>::print(int level);
template void qsTensorTrain<std::complex<double>, 1>::print(int level);
template void qsTensorTrain<std::complex<double>, 2>::print(int level);
//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
// Basic arithmetic operations
template <typename T, unsigned N>
qsTensorTrain<T, N>& qsTensorTrain<T, N>::operator = (const qsTensorTrain<T, N>& other){
  if(this!=&other && other.tensors_allocated){
    _id = other._id;
    freeTensors();
    setLength(other.length);
    setPhysicalDim(other.phy_dim);
    bond_dims = other.bond_dims;
    totalQ = other.totalQ;
    phy_qn = other.phy_qn;
    center = other.center;
    A.clear();
    A = other.A;
    tensors_allocated = true;
  }
  return *this;
}
template qsTensorTrain<double, 1>& qsTensorTrain<double, 1>::operator = (const qsTensorTrain<double, 1>& other);
template qsTensorTrain<double, 2>& qsTensorTrain<double, 2>::operator = (const qsTensorTrain<double, 2>& other);
template qsTensorTrain<std::complex<double>, 1>& qsTensorTrain<std::complex<double>, 1>::operator = (const qsTensorTrain<std::complex<double>, 1>& other);
template qsTensorTrain<std::complex<double>, 2>& qsTensorTrain<std::complex<double>, 2>::operator = (const qsTensorTrain<std::complex<double>, 2>& other);


template <typename T, unsigned N>
qsTensorTrain<T, N>& qsTensorTrain<T, N>::operator = (qsTensorTrain<T, N>&& other){
  if(this!=&other && other.tensors_allocated){
    _id = other._id;
    freeTensors();
    setLength(other.length);
    setPhysicalDim(other.phy_dim);
    bond_dims = std::move(other.bond_dims);
    totalQ = other.totalQ;
    phy_qn = std::move(other.phy_qn);
    center = other.center;
    A.clear();
    A = std::move(other.A);
    tensors_allocated = true;
    other.tensors_allocated = false;
  }
  return *this;
}
template qsTensorTrain<double, 1>& qsTensorTrain<double, 1>::operator = (qsTensorTrain<double, 1>&& other);
template qsTensorTrain<double, 2>& qsTensorTrain<double, 2>::operator = (qsTensorTrain<double, 2>&& other);
template qsTensorTrain<std::complex<double>, 1>& qsTensorTrain<std::complex<double>, 1>::operator = (qsTensorTrain<std::complex<double>, 1>&& other);
template qsTensorTrain<std::complex<double>, 2>& qsTensorTrain<std::complex<double>, 2>::operator = (qsTensorTrain<std::complex<double>, 2>&& other);


template <typename T, unsigned N>
qsTensorTrain<T, N>& qsTensorTrain<T, N>::operator = (const qTensorTrain<T, N>& other){ //convert
  if(other.tensors_allocated){
    _id = other._id;
    freeTensors();
    setLength(other.length);
    setPhysicalDim(other.phy_dim);
    bond_dims = other.bond_dims;
    totalQ = other.totalQ;
    phy_qn = other.phy_qn;
    center = other.center;
    A.clear();
    A.resize(other.A.size());
    for(int i=0;i<A.size();i++)
      A[i] = other.A[i];
    tensors_allocated = true;
  }
  return *this;
}
template qsTensorTrain<double, 1>& qsTensorTrain<double, 1>::operator = (const qTensorTrain<double, 1>& other);
template qsTensorTrain<double, 2>& qsTensorTrain<double, 2>::operator = (const qTensorTrain<double, 2>& other);
template qsTensorTrain<std::complex<double>, 1>& qsTensorTrain<std::complex<double>, 1>::operator = (const qTensorTrain<std::complex<double>, 1>& other);
template qsTensorTrain<std::complex<double>, 2>& qsTensorTrain<std::complex<double>, 2>::operator = (const qTensorTrain<std::complex<double>, 2>& other);


template <typename T, unsigned N>
qsTensorTrain<T, N>& qsTensorTrain<T, N>::operator = (qTensorTrain<T, N>&& other){
  if(other.tensors_allocated){
    _id = other._id;
    freeTensors();
    setLength(other.length);
    setPhysicalDim(other.phy_dim);
    bond_dims = std::move(other.bond_dims);
    totalQ = other.totalQ;
    phy_qn = std::move(other.phy_qn);
    center = other.center;
    A.clear();
    //A = std::move(other.A);
    A.resize(other.A.size());
    for(int i=0;i<A.size();i++)
      A[i] = other.A[i];
    tensors_allocated = true;
    other.tensors_allocated = false;
  }
  return *this;
}
template qsTensorTrain<double, 1>& qsTensorTrain<double, 1>::operator = (qTensorTrain<double, 1>&& other);
template qsTensorTrain<double, 2>& qsTensorTrain<double, 2>::operator = (qTensorTrain<double, 2>&& other);
template qsTensorTrain<std::complex<double>, 1>& qsTensorTrain<std::complex<double>, 1>::operator = (qTensorTrain<std::complex<double>, 1>&& other);
template qsTensorTrain<std::complex<double>, 2>& qsTensorTrain<std::complex<double>, 2>::operator = (qTensorTrain<std::complex<double>, 2>&& other);

template <typename T, unsigned N>
qsTensorTrain<T, N>& qsTensorTrain<T, N>::operator *= (const T c){
  assert(tensors_allocated);
  A[0] *= c;
  return *this;
}
template qsTensorTrain<double, 1>& qsTensorTrain<double, 1>::operator *= (const double c);
template qsTensorTrain<double, 2>& qsTensorTrain<double, 2>::operator *= (const double c);
template qsTensorTrain<std::complex<double>, 1>& qsTensorTrain<std::complex<double>, 1>::operator *= (const  std::complex<double> c);
template qsTensorTrain<std::complex<double>, 2>& qsTensorTrain<std::complex<double>, 2>::operator *= (const  std::complex<double> c);


template <typename T, unsigned N>
qsTensorTrain<T, N>& qsTensorTrain<T, N>::operator /= (const T c){
  assert(tensors_allocated);
  A[0] /= c;
  return *this;
}
template qsTensorTrain<double, 1>& qsTensorTrain<double, 1>::operator /= (const double c);
template qsTensorTrain<double, 2>& qsTensorTrain<double, 2>::operator /= (const double c);
template qsTensorTrain< std::complex<double>, 1>& qsTensorTrain<std::complex<double>, 1>::operator /= (const  std::complex<double> c);
template qsTensorTrain< std::complex<double>, 2>& qsTensorTrain<std::complex<double>, 2>::operator /= (const  std::complex<double> c);


template <typename T, unsigned N>
qsTensorTrain<T, N> qsTensorTrain<T, N>::operator * (const T c) const{
  assert(tensors_allocated);
  qsTensorTrain<T, N> t;
  t = *this;
  t.A[0] *= c;
  return t;
}
template qsTensorTrain<double, 1> qsTensorTrain<double, 1>::operator * (const double c) const;
template qsTensorTrain<double, 2> qsTensorTrain<double, 2>::operator * (const double c) const;
template qsTensorTrain< std::complex<double>, 1> qsTensorTrain<std::complex<double>, 1>::operator * (const  std::complex<double> c) const;
template qsTensorTrain< std::complex<double>, 2 > qsTensorTrain<std::complex<double>, 2>::operator * (const  std::complex<double> c) const;


template <typename T, unsigned N>
qsTensorTrain<T, N> qsTensorTrain<T, N>::operator / (const T c) const{
  assert(tensors_allocated);
  qsTensorTrain<T, N> t;
  t = *this;
  t.A[0] /= c;
  return t;
}
template qsTensorTrain<double, 1> qsTensorTrain<double, 1>::operator / (const double c) const;
template qsTensorTrain<double, 2> qsTensorTrain<double, 2>::operator / (const double c) const;
template qsTensorTrain< std::complex<double>, 1> qsTensorTrain<std::complex<double>, 1>::operator / (const  std::complex<double> c) const;
template qsTensorTrain< std::complex<double>, 2> qsTensorTrain<std::complex<double>, 2>::operator / (const  std::complex<double> c) const;


template <typename T, unsigned N>
void qsTensorTrain<T, N>::save(std::string fn){
  assert(tensors_allocated);
  ezh5::File fh5 (fn, H5F_ACC_TRUNC);
  fh5["length"] = length;
  fh5["phy_dim"] = phy_dim;
  fh5["bond_dims"] = bond_dims;
  fh5["center"] = center;
  fh5["totalQ"] = totalQ;
  fh5["phy_qn"] = phy_qn;
  fh5["id"] = _id;
  for (size_t i = 0; i < length; i++){
    ezh5::Node nd = fh5["Tensor"+to_string(i)];
    A[i].save(nd);
  }
}
template void qsTensorTrain<double, 1>::save(std::string fn);
template void qsTensorTrain<double, 2>::save(std::string fn);
template void qsTensorTrain<std::complex<double>, 1>::save(std::string fn);
template void qsTensorTrain<std::complex<double>, 2>::save(std::string fn);


template <typename T, unsigned N>
void qsTensorTrain<T, N>::load(std::string fn){
  freeTensors();
  ezh5::File fh5 (fn, H5F_ACC_RDONLY);
  fh5["length"] >> length;
  fh5["phy_dim"] >> phy_dim;
  fh5["bond_dims"] >> bond_dims;
  fh5["center"] >> center;
  fh5["totalQ"] >> totalQ;
  fh5["phy_qn"] >> phy_qn;
  fh5["id"] >> _id;
  allocateTensors();
  for (size_t i = 0; i < length; i++){
    ezh5::Node nd = fh5["Tensor"+to_string(i)];
    A[i].load(nd);
  }
}
template void qsTensorTrain<double, 1>::load(std::string fn);
template void qsTensorTrain<double, 2>::load(std::string fn);
template void qsTensorTrain<std::complex<double>, 1>::load(std::string fn);
template void qsTensorTrain<std::complex<double>, 2>::load(std::string fn);


template <typename T, unsigned N>
void qsTensorTrain<T, N>::save(ezh5::Node& fh5){
  assert(tensors_allocated);
  fh5["length"] = length;
  fh5["phy_dim"] = phy_dim;
  fh5["bond_dims"] = bond_dims;
  fh5["center"] = center;
  fh5["totalQ"] = totalQ;
  fh5["phy_qn"] = phy_qn;
  fh5["id"] = _id;
  for (size_t i = 0; i < length; i++){
    ezh5::Node nd = fh5["Tensor"+to_string(i)];
    A[i].save(nd);
  }
}
template void qsTensorTrain<double, 1>::save(ezh5::Node& fh5);
template void qsTensorTrain<double, 2>::save(ezh5::Node& fh5);
template void qsTensorTrain<std::complex<double>, 1>::save(ezh5::Node& fh5);
template void qsTensorTrain<std::complex<double>, 2>::save(ezh5::Node& fh5);


template <typename T, unsigned N>
void qsTensorTrain<T, N>::load(ezh5::Node& fh5){
  freeTensors();
  fh5["length"] >> length;
  fh5["phy_dim"] >> phy_dim;
  fh5["bond_dims"] >> bond_dims;
  fh5["center"] >> center;
  fh5["totalQ"] >> totalQ;
  fh5["phy_qn"] >> phy_qn;
  fh5["id"] >> _id;
  allocateTensors();
  for (size_t i = 0; i < length; i++){
    ezh5::Node nd = fh5["Tensor"+to_string(i)];
    A[i].load(nd);
  }
}
template void qsTensorTrain<double, 1>::load(ezh5::Node& fh5);
template void qsTensorTrain<double, 2>::load(ezh5::Node& fh5);
template void qsTensorTrain<std::complex<double>, 1>::load(ezh5::Node& fh5);
template void qsTensorTrain<std::complex<double>, 2>::load(ezh5::Node& fh5);


template <typename T, unsigned N>
void qsTensorTrain<T, N>::rc(){
  assert(tensors_allocated);
  qstensor<T> U,V,S;
  for (size_t i = length-1; i > 0; i--) {
    vector<qtensor_index> left;
    vector<qtensor_index> right;
    qtensor_index mid;
    string LinkName = "ID"+to_string(_id)+"Link"+to_string(i);
    // Separate qtensor_index
    for (size_t j = 0; j < A[i].rank; j++) {
      string idx_name = A[i].idx_set[j].name();
      if (idx_name == LinkName) {
        left.push_back(A[i].idx_set[j]);
        mid = A[i].idx_set[j];
      }else{
        right.push_back(A[i].idx_set[j]);
      }
    }
    // SVD
    svd(A[i],left,right,U,V,S,MoveFromRight);
    unsigned S_size = 0;
    for(auto & block: S._block) S_size+= block.get_tot_size(false);
    bond_dims[i] = S_size;
    A[i] = V;
    A[i].idx_set[0].rename(mid.name());
    A[i].idx_set[0].prime(mid.level()-A[i].idx_set[0].level());
    A[i-1] = std::move(A[i-1]*U);
    A[i-1].idx_set.back().rename(mid.name());
    A[i-1].idx_set.back().prime(mid.level()-A[i-1].idx_set.back().level());
  }
  center = 0;
}
template void qsTensorTrain<double, 1>::rc();
template void qsTensorTrain<double, 2>::rc();
template void qsTensorTrain<std::complex<double>, 1>::rc();
template void qsTensorTrain<std::complex<double>, 2>::rc();


template <typename T, unsigned N>
void qsTensorTrain<T, N>::lc(){
  assert(tensors_allocated);
  qstensor<T> U,V,S;
  for (size_t i = 0; i < length-1; i++) {
    vector<qtensor_index> left;
    vector<qtensor_index> right;
    qtensor_index mid;
    string LinkName = "ID"+to_string(_id)+"Link"+to_string(i+1);
    // Separate qtensor_index
    for (size_t j = 0; j < A[i].rank; j++) {
      string idx_name = A[i].idx_set[j].name();
      if (idx_name == LinkName) {
        right.push_back(A[i].idx_set[j]);
        mid = A[i].idx_set[j];
      }else{
        left.push_back(A[i].idx_set[j]);
      }
    }
    // SVD
    svd(A[i],left,right,U,V,S,MoveFromLeft);
    unsigned S_size = 0;
    for(auto & block: S._block) S_size+= block.get_tot_size(false);
    bond_dims[i+1] = S_size;
    A[i] = U;
    A[i].idx_set.back().rename(mid.name());
    A[i].idx_set.back().prime(mid.level()-A[i].idx_set.back().level());
    A[i+1] = std::move(V*A[i+1]);
    A[i+1].idx_set[0].rename(mid.name());
    A[i+1].idx_set[0].prime(mid.level()-A[i+1].idx_set[0].level());
  }
  center = length-1;
}
template void qsTensorTrain<double, 1>::lc();
template void qsTensorTrain<double, 2>::lc();
template void qsTensorTrain<std::complex<double>, 1>::lc();
template void qsTensorTrain<std::complex<double>, 2>::lc();


template <typename T, unsigned N>
void qsTensorTrain<T, N>::normalize(){
  assert(tensors_allocated);
  if(center == -1){
    double nm = norm();
    A[0] /= nm;
  } else{
    double nm = A[center].norm();
    A[center] /= nm;
  }
}
template void qsTensorTrain<double, 1>::normalize();
template void qsTensorTrain<double, 2>::normalize();
template void qsTensorTrain<std::complex<double>, 1>::normalize();
template void qsTensorTrain<std::complex<double>, 2>::normalize();


template <typename T, unsigned N>
double qsTensorTrain<T, N>::norm(){
  assert(tensors_allocated);
  qsTensorTrain<T, N> t;
  t = *this;
  t.rc();
  return t.A[0].norm();
}
template double qsTensorTrain<double, 1>::norm();
template double qsTensorTrain<double, 2>::norm();
template double qsTensorTrain<std::complex<double>, 1>::norm();
template double qsTensorTrain<std::complex<double>, 2>::norm();


//---------------------------------------------------------------------------
// Adjust canonical center
template <typename T, unsigned N>
double qsTensorTrain<T, N>::position(int site){
  assert(tensors_allocated);
  assert(site>=0 && site<int(length));
  // Initialize center position if not in canonical form
  if(center == -1) rc();
  qstensor<T> U,V,S;
  while( center!=site ){
    vector<qtensor_index> left;
    vector<qtensor_index> right;
    qtensor_index mid;
    if(center>site){
      // Move center to the left
      string LinkName = "ID"+to_string(_id)+"Link"+to_string(center);
      for (size_t j = 0; j < A[center].rank; j++) {
        string idx_name = A[center].idx_set[j].name();
        if (idx_name == LinkName) {
          left.push_back(A[center].idx_set[j]);
          mid = A[center].idx_set[j];
        }else{
          right.push_back(A[center].idx_set[j]);
        }
      }
      svd(A[center],left,right,U,V,S,MoveFromRight);
      unsigned S_size = 0;
      for(auto & block: S._block) S_size+= block.get_tot_size(false);
      bond_dims[center] = S_size;
      A[center] = V;
      A[center].idx_set[0].rename(mid.name());
      A[center].idx_set[0].prime(mid.level()-A[center].idx_set[0].level());
      A[center-1] = std::move(A[center-1]*U);
      A[center-1].idx_set.back().rename(mid.name());
      A[center-1].idx_set.back().prime(mid.level()-A[center-1].idx_set.back().level());
      --center;
    }else if(center<site){
      // Move center to the right
      string LinkName = "ID"+to_string(_id)+"Link"+to_string(center+1);
      for (size_t j = 0; j < A[center].rank; j++) {
        string idx_name = A[center].idx_set[j].name();
        if (idx_name == LinkName) {
          right.push_back(A[center].idx_set[j]);
          mid = A[center].idx_set[j];
        }else{
          left.push_back(A[center].idx_set[j]);
        }
      }
      svd(A[center],left,right,U,V,S,MoveFromLeft);
      unsigned S_size = 0;
      for(auto & block: S._block) S_size+= block.get_tot_size(false);
      bond_dims[center+1] = S_size;
      A[center] = U;
      A[center].idx_set.back().rename(mid.name());
      A[center].idx_set.back().prime(mid.level()-A[center].idx_set.back().level());
      A[center+1] = std::move(V*A[center+1]);
      A[center+1].idx_set[0].rename(mid.name());
      A[center+1].idx_set[0].prime(mid.level()-A[center+1].idx_set[0].level());
      ++center;
    }
  }
  if(S._initted == false) return 0.; //didn't actually move
  double vNEE = calcEntropy(S);
  return vNEE;
}
template double qsTensorTrain<double, 1>::position(int site);
template double qsTensorTrain<double, 2>::position(int site);
template double qsTensorTrain<std::complex<double>, 1>::position(int site);
template double qsTensorTrain<std::complex<double>, 2>::position(int site);


#endif
