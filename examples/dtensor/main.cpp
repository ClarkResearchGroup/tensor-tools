#include "../../dtensor/tensor_index.cpp"
#include "../../dtensor/dtensor.cpp"
#include "../../util/ezh5.cpp"
#include "../../util/timer.h"
int main(int argc, char const *argv[]) {
  //------------------------------------------
  std::cout << "Hello Tensor!" << '\n';
  // dtensor<double> A({3,4,5});
  // std::cout << A._T.stride(0) << '\n';
  // std::cout << A._T.stride(1) << '\n';
  // std::cout << A._T.stride(2) << '\n';
  // std::cout << A.idx_set[0]._name << '\n';
  // std::cout << A.idx_set[1]._name << '\n';
  // std::cout << A.idx_set[2]._name << '\n' << '\n';
  // std::cout << (A.idx_set[0] == A.idx_set[1]) << '\n';
  //
  // dtensor<double> B = std::move(A);
  // std::cout << B._T.stride(0) << '\n';
  // std::cout << B._T.stride(1) << '\n';
  // std::cout << B._T.stride(2) << '\n';
  // std::cout << B.idx_set[0]._name << '\n';
  // std::cout << B.idx_set[1]._name << '\n';
  // std::cout << B.idx_set[2]._name << '\n';
  // std::cout << B.idx_set[0]._size << '\n';
  // std::cout << B.idx_set[1]._size << '\n';
  // std::cout << B.idx_set[2]._size << '\n';
  // std::cout << B.idx_set[0]._type << '\n';
  // std::cout << B.idx_set[1]._type << '\n';
  // std::cout << B.idx_set[2]._type << '\n';

  // tensor<double> T1({4,5});
  // tensor<double> T2({5,4});
  // tensor<double> T3;
  // Mxd mT1(4,5), mT2(5,4);
  // mT1.setRandom(); mT2.setRandom();
  // std::copy(mT1.data(),mT1.data()+mT1.size(),T1.data());
  // std::copy(mT2.data(),mT2.data()+mT2.size(),T2.data());
  // std::vector<tblis::label_type> idx_A(2);
  // idx_A[0] = 'a';
  // idx_A[1] = 'b';
  // std::vector<tblis::label_type> idx_B(2);
  // idx_B[0] = 'b';
  // idx_B[1] = 'a';
  // double h = 0.0;
  // tblis::dot(T1,idx_A.data(),T2,idx_B.data(),h);
  // std::cout << h << '\n';
  // std::cout << (mT1*mT2).trace() << '\n' << '\n';
  // for (int i = 0; i < 20; i++) {
  //   std::cout << T1.data()[i] << " ";
  // }
  // std::cout << '\n';
  // tblis::scale(2.0, T1, idx_A.data());
  // for (int i = 0; i < 20; i++) {
  //   std::cout << T1.data()[i] << " ";
  // }
  // std::cout << '\n';
  //
  // T3.reset(std::move(T1));
  // for (int i = 0; i < 20; i++) {
  //   std::cout << T3.data()[i] << " ";
  // }
  // std::cout << '\n';

  // tensor<double> T1({4,5},0.1);
  // tensor<double> T2({6,7},0.2);
  // tensor<double> T3({4,5,6,7});
  // std::vector<tblis::label_type> idx_A(2);
  // idx_A[0] = 'a';
  // idx_A[1] = 'b';
  // std::vector<tblis::label_type> idx_B(2);
  // idx_B[0] = 'c';
  // idx_B[1] = 'd';
  // std::vector<tblis::label_type> idx_C(4);
  // idx_C[0] = 'a';
  // idx_C[1] = 'b';
  // idx_C[2] = 'c';
  // idx_C[3] = 'd';
  // tblis::mult(1.0,T1,idx_A.data(),T2,idx_B.data(),1.0,T3,idx_C.data());
  // for (int i = 0; i < 20*6*7; i++) {
  //   std::cout << T3.data()[i] << " ";
  // }
  // std::cout << '\n';

  // T3.reset(std::move(T1));
  // for (int i = 0; i < 20; i++) {
  //   std::cout << T3.data()[i] << " ";
  // }
  // std::cout << '\n';

  // dtensor<double> A({3,4,5},{"a","b","c"}); A.setRandom(); A.print();
  // dtensor<double> B({4,7,5},{"b","d","c"}); B.setRandom(); B.print();
  // dtensor<double> C = (B*A); C.print(1);
  // dtensor<double> D({3,7},{"a","d"}); D.setRandom(); D.print(1);
  // double res = C.contract(D);
  // std::cout << res << '\n';

  // dtensor<double> A({4,4},{"a","a"},{Site,Site},{0,1});
  // A.setRandom(); A.print(1);
  // dtensor<double> B = A.diagonal(); B.print(1);
  // A.save("test.h5");

  dtensor<double> A({2,3});
  A.load("test.h5");
  A.mapPrime(0,2);
  A.normalize();
  A.print(1);
  //------------------------------------------
  return 0;
}
