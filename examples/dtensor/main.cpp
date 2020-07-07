/*
 * Copyright 2020 Xiongjie Yu, Ryan Levy, and Bryan K. Clark
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

#include "../../dtensor/dtensor_index.cpp"
#include "../../dtensor/dtensor_index_op.cpp"
#include "../../dtensor/dtensor.cpp"
#include "../../dtensor/dtensor_view.cpp"
#include "../../dtensor/dtensor_op.cpp"
#include "../../dtensor/big_dtensor.cpp"
#include "../../linalg/lapack_wrapper.cpp"
#include "../../linalg/tensor_cg.cpp"
#include "../../linalg/tensor_davidson.cpp"
#include "../../util/ezh5.cpp"
#include "../../util/timer.h"
int main(int argc, char const *argv[]) {
  //------------------------------------------
  std::cout << "Hello Tensor!" << '\n';

  unsigned s1=100, s2=100;
  Eigen::MatrixXd Id(s1,s2); Id.setIdentity();
  Eigen::MatrixXd mA(s1,s2); mA.setRandom();
  Eigen::MatrixXd tp = mA + mA.transpose();
  mA = tp;
  Eigen::MatrixXd mB(s1,1);   mB.setRandom();
  // std::cout << mA << '\n' << '\n';
  // std::cout << mB.transpose() << '\n' << '\n';
  std::cout << (mA*mB).transpose() << '\n' << '\n';
  std::cout << mB.transpose()*mA*mB << '\n' << '\n';
  std::cout << mA.fullPivLu().solve(mB).transpose() << '\n' << '\n';
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(mA);
  std::cout<<es.eigenvalues().transpose() << '\n' << '\n';
  std::cout<<es.eigenvectors().col(0).transpose() << '\n' << '\n';
  std::cout<<es.eigenvectors().col(s1-1).transpose() << '\n' << '\n';

  dtensor<double> A({s1,s2},{"a","a"},{Link,Link},{1,0},mA.data());
  dtensor<double> B({s1},{"a"},{Link},{0},mB.data());
  A.print();
  B.print();
  big_dtensor<double> M;
  M.addMid(&A);
  dtensor<double> C = M.product(B);
  C.print(1);
  std::cout << M.expec(B) << '\n';
  C.setRandom();
  tensor_CG(M, C, B, 300, 1e-10, false);
  C.print(1);
  std::cout << tensor_davidson(M, B, 5, 100, 1e-10, 'S') << '\n';
  B.print(1);
  //------------------------------------------
  return 0;
}
