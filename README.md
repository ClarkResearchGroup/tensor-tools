<img src="https://github.com/ClarkResearchGroup/tensor-tools/raw/parallel/logo2.png" width="250px" align="right" alt="tensor-tools logo" />  

# tensor-tools parallel 


### What is it?

The parallel verison of `tensor-tools` is a parallel tensor-network library currently customized for performing DMRG (with quantum numbers) in parallel. The parallelization is done via MPI over nodes, using [CTF tensors](https://github.com/cyclops-community/ctf/). It is particularly powerful for spin-systems or those without conserved U(1) symmetries.  
Not only can it reach bond-dimensions which are inaccessible serially, but (if we extrapolate what a serial code could do on a big-enough machine) it achieves a 10x speedup at no extra cost in node-hours.  See [arxiv:2007.05540](https://arxiv.org/abs/2007.05540).  

### Authors

The parallel version of tensor-tools was developed by [Ryan Levy](https://ryanlevy.github.io/) under the guidance of [Bryan Clark](http://clark.physics.illinois.edu/) with CTF and other help from [Edgar Solomonik](http://solomonik.cs.illinois.edu/).  It is built on top of the [serial version of tensor-tools](https://github.com/ClarkResearchGroup/tensor-tools/tree/serial) which was written (circa 2017) by Xiongjie Yu (under the guidance of Clark); the serial code now also has additional contributions/bug-fixes from Ryan Levy.


## Installation and Running

Requirements:
- [CTF](https://github.com/cyclops-community/ctf/)
  - compiled with lapack/scalapack support; hptt support preferred 
  - example configure line: `./configure --no-dynamic --build-hptt --with-scalapack' 'LIB_PATH=-L${MKLROOT}/lib/intel64/' 'LIBS=-lmkl_scalapack_lp64 -Wl,--start-group -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_blacs_intelmpi_lp64 -Wl,--end-group -lpthread -lm'`
  - If using sparse tensors, intel MKL is strongly suggested
- HDF5 
- [Eigen3](http://eigen.tuxfamily.org/index.php?title=Main_Page)
- C++11 support (GCC 6.3.0/8.2.0, Intel 18.0.2/18.0.3.222 compilers tested)

After installing the required libraries, simply make a copy of `Makefile.in` and adjust the variables to their pointed location. Examples of a linux system, Blue Waters, and Stampede2 are provided. Then build (`make -j`) which will create an executable `run` based on the `project/main.cpp` file. 

To execute, launch with MPI using `mpirun -n [# of procs] run [inputFile]` and follow CTF good practice such as `export OMP_NUM_THREADS=1` and `CTF_PPN=[# of procs/node]`.  

## Conversion

Conversion to and from ITensor files is fully supported for v2 and v3.1.x (to be determined). See [link repository here]. 

-----------------

### External Code
We utilize both the AutoMPO system from [ITensor](https://github.com/ITensor/ITensor/) and some of their lattice code, to have better interoperability between the two. See [EXTERNAL.md](EXTERNAL.md) for more.  

Some of the general tensor-tools interface/API was also inspired by ITensor.  
