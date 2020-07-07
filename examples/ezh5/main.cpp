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


#include <iostream>
#include <vector>
#include <string>
#include <Eigen/Dense>

#include "../../util/ezh5.cpp"

typedef Eigen::MatrixXd Mxd;

using namespace std;

void save_by_handle(ezh5::Node& f, double val){
  f["pi"] = val;
}

void read_by_handle(ezh5::Node& f, double& val){
  f["pi"] >> val;
}

int main(int argc, char const *argv[]) {
  cout<<"//------------------------------------"<<endl;
  cout<<"This program tests the simplified HDF5 interface provided by ezh5."<<endl;
  cout<<"It will write/read a random vector and a random Eigen3 matrix to/from files."<<endl;
  cout<<"//------------------------------------"<<endl;
  //------------------------------------
  int reading_mode = 1;
  int writing_mode = 0;
  int mode = 0;
  //------------------------------------
  int N;
  cout<<"Please specify write (0) or read (1): ";
  cin>>mode;
  //------------------------------------
  string fn1 = "test_stl_vector.h5";
  string fn2 = "test_Eigen3_matrix.h5";
  //------------------------------------
  if (mode == writing_mode) {
    cout<<"Please specify vector/matrix size: ";
    cin>>N;
    //------------------------------------
    cout<<"//------------------------------------"<<endl;
    vector<double> v(N);
    for (int i = 0; i < N; i++) v[i] = drand48();
    cout<<"Random vector of length "<<N<<":"<<endl;
    for(auto val : v) cout<<val<<" ";
    cout<<endl;
    ezh5::File f1 (fn1, H5F_ACC_TRUNC);
    f1["my_random_vector"] = v;
    cout<<"Vector is written to '"<<fn1<<"'."<<endl;
    cout<<"If hdf5 tools are installed, use 'h5dump "<<fn1<<"' to display the file content in terminal."<<endl;
    double pi = 3.141592653;
    ezh5::Node f1n = f1["const"];
    save_by_handle(f1n, pi);
    //------------------------------------
    cout<<"//------------------------------------"<<endl;
    Mxd M(N,N); M.setRandom();
    cout<<"Random Eigen3 matrix of size "<<N<<"x"<<N<<":"<<endl;
    cout<<M<<endl;
    ezh5::File f2 (fn2, H5F_ACC_TRUNC);
    f2["my_random_matrix"] = M;
    cout<<"Matrix is written to '"<<fn2<<"'."<<endl;
    cout<<"If hdf5 tools are installed, use 'h5dump "<<fn2<<"' to display the file content in terminal."<<endl;
    cout<<"//------------------------------------"<<endl;
  }else if (mode == reading_mode) {
    //------------------------------------
    cout<<"//------------------------------------"<<endl;
    vector<double> v;
    ezh5::File f1 (fn1, H5F_ACC_RDONLY);
    f1["my_random_vector"] >> v;
    for(auto val : v) cout<<val<<" ";
    cout<<endl;
    cout<<"Vector is read from '"<<fn1<<"'."<<endl;
    double pi = 0;
    ezh5::Node f1n = f1["const"];
    read_by_handle(f1n, pi);
    std::cout << "Pi = " << pi << '\n';
    //------------------------------------
    cout<<"//------------------------------------"<<endl;
    Mxd M;
    ezh5::File f2 (fn2, H5F_ACC_RDONLY);
    f2["my_random_matrix"] >> M;
    cout<<M<<endl;
    cout<<"Matrix is read from '"<<fn2<<"'."<<endl;
    cout<<"//------------------------------------"<<endl;
  }
  //------------------------------------
  return 0;
}
