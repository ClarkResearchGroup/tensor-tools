# tensor-tools

-----------------
Dependency
```C++
//------------------------------------------------------------------------------
// External library 1 -- TBLIS -- necessary
//------------------------------------------------------------------------------
// By Devin Matthews, et al.
// https://github.com/devinamatthews/tblis
//          **Make sure to use v1.1.2**
// Intro:
// TBLIS is a library and framework for performing tensor operations,
// especially tensor contraction, using efficient native algorithms.
// In other words, TBLIS does not need tensor transposition to perform tensor
// contraction.
// In comparison, a popular but not optimal tensor contraction strategy consists
// of 3 steps: transposition, d/s/zgemm, transposition.
// However, TBLIS does not support sparse tensor (at least for now).
#include "tblis.h"
typedef vector<tblis::len_type>    len_vec;
typedef vector<tblis::label_type>  lab_vec;
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// External library 2 -- Eigen -- needed, but only for conducting easier tests
//------------------------------------------------------------------------------
// http://eigen.tuxfamily.org
#include <Eigen/Dense>
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// External library 3 -- HDF5 -- required when saving/loading data
//------------------------------------------------------------------------------
// Check https://support.hdfgroup.org/HDF5/ for instructions on
// how to install HDF5
#include "hdf5.h"
//------------------------------------------------------------------------------
// For easier user interface of HDF5 functions, a wrapper called "EZH5" is used.
// By M. Chen
// https://github.com/mileschen360/ezh5
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// External library 4 -- HPTT -- optional
//------------------------------------------------------------------------------
// By Paul Springer
// https://github.com/springer13/hptt
// Intro:
// HPTT is a high-performance C++ library for out-of-place tensor transpositions.
// Notice:
// Tensor transposition is used to reshape tensors before sending them to QR/SVD.
// If enabled, tensor transposition calls will be sent to HPTT driver.
// If not enabled, tensor transposition will be sent to a not optimal function.
#ifdef USE_HPTT
#include <hptt.h>
#endif
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// Single Abelian quantum number
// first: qn;    second: dimension
//------------------------------------------------------------------------------
typedef pair<int, unsigned> quantum_number;
//------------------------------------------------------------------------------
```

-----------------
I. Dense tensor interfaces
```C++
// Create a dtensor<T> object by providing sizes, names, types (Link, Site), prime level (unsigned)
dtensor<double> A({4,8,6}, {"a","b","c"}, {Link,Link}, {0,0});
A.setRandom();
// normalize
A.normalize();
// print has different levels of information available
A.print(1);
// Create a dtensor<T> object by providing sizes, names (default to Link type and 0 prime level)
dtensor<double> B({6,7,4}, {"c","d","a"});
B.setRandom();
B.print(1);
// Repeated indices (same size, name, type and prime level) are contracted in multiplication
// dtensor<T> supports move semantics, which removes unnecessary copy operation
dtensor<double> C = std::move(A*B);
C.print(1);
```

-----------------
II. Quantum numbered tensor interfaces
```C++
// Create a qtensor<T> object by providing arrows, names (default to Link type and 0 prime level)
qtensor< std::complex<double> > A({Inward,Inward,Outward} ,{"a", "s1", "b"});
// Fill in the information on the quantum numbers (qn, qdim)
A.addQNtoIndex(0, std::make_pair(1, 4));
A.addQNtoIndex(0, std::make_pair(2, 5));
A.addQNtoIndex(1, std::make_pair(0, 1));
A.addQNtoIndex(1, std::make_pair(1, 1));
A.addQNtoIndex(2, std::make_pair(1, 6));
A.addQNtoIndex(2, std::make_pair(2, 7));
A.addQNtoIndex(2, std::make_pair(3, 8));
// Create blocks that have zero total quantum number
A.initBlock();
A.setRandom();
// print has different levels of information available
A.print(1);

qtensor< std::complex<double> > B({Inward,Inward,Outward} ,{"b", "s2", "c"});
B.addQNtoIndex(0, std::make_pair(1, 6));
B.addQNtoIndex(0, std::make_pair(2, 7));
B.addQNtoIndex(0, std::make_pair(3, 8));
B.addQNtoIndex(1, std::make_pair(0, 1));
B.addQNtoIndex(1, std::make_pair(1, 1));
B.addQNtoIndex(2, std::make_pair(1, 7));
B.addQNtoIndex(2, std::make_pair(2, 8));
B.addQNtoIndex(2, std::make_pair(3, 9));
B.addQNtoIndex(2, std::make_pair(4, 10));
B.initBlock();
B.setRandom();
B.print(1);

// qtensor<T> multiplication
qtensor< std::complex<double> > C = std::move(A*B);
C.print(1);
```

-----------------
III. DMRG interfaces
```C++
unsigned N = 30;
// Set sites
spinhalf sites(N);
// Create product string
vector<string> ps;
for (size_t i = 0; i < sites.N(); i++) {
  if(i%2==0)
    ps.push_back("Dn");
  else
    ps.push_back("Up");
}

{
  std::cout << "DMRG without quantum number." << '\n';
  MPS<double> psi(&sites);
  psi.setRandom();
  psi.normalize();

  MPO<double> H;
  Heisenberg<double> HB(&sites);
  HB.buildHam(H);

  int nsweeps = 20;
  int maxm = 60;
  double cutoff = 1e-8;
  dmrg(psi, H, nsweeps, maxm, cutoff);

  psi.print(0);
  std::cout << psiHphi(psi,H,psi) << '\n';
}

{
  std::cout << "DMRG with total Sz quantum number." << '\n';
  qMPS<double> psi(&sites,ps);
  psi.normalize();

  qMPO<double> H;
  Heisenberg<double> HB(&sites);
  HB.buildHam(H);

  int nsweeps = 20;
  int maxm = 60;
  double cutoff = 2e-11;
  dmrg(psi, H, nsweeps, maxm, cutoff);

  psi.print(0);
  std::cout << psiHphi(psi,H,psi) << '\n';
}
```
