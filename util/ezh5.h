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

#ifndef EZH5_H
#define EZH5_H

#include <iostream>
#include "hdf5.h"
#include <Eigen/Dense>
#include <vector>

namespace ezh5{
  ///////////////////////////
  class ID{
    public:
    hid_t id;
    ID(): id(-1) {}
    ID(hid_t id_in) : id(id_in) {}
    ~ID(){}
  };
  ///////////////////////////
  class Node :public ID
  {
  public:
    hid_t pid;
    std::string path;
    bool verbose_;

    Node(){verbose_ = false;}
    Node(hid_t pid_in, const std::string& path_in, bool verbose=false): ID(-1), pid(pid_in), path(path_in){
    	verbose_ = verbose;
    	if(verbose_) std::cout<<"creating "<<path<<std::endl;
    }
    ////////
    Node& operator()(const std::string& path){
    	return *this;
    }
    ////////
    Node operator[] (const std::string& path_more);
    ////////
    template<typename T> Node& operator =  (T val);
    template<typename T> Node& operator >> (T& val);
    ////////
    template<typename T> Node& operator =  (std::vector<T>& vec);
    template<typename T> Node& operator >> (std::vector<T>& vec);
    ////////
    template<typename T> Node& operator =  (Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& m);
    template<typename T> Node& operator >> (Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& m);
    ////////
    ~Node(){
    	if(verbose_) std::cout<<"closing "<<path<<std::endl;
    	if(this->id>0)
    	{
    		H5Gclose(this->id);
    		this->id = -1;
    	}
    }
  };
  ///////////////////////////
  class File : public Node
  {
  public:
    File(const std::string& path, unsigned flags, bool verbose = false);

    ~File()
    {
    	H5Fclose(this->id);
    	this->id = -1;
    }
  };
  ///////////////////////////
}


#endif
