# CC = icpc
# CC_FLAGS = -O3 -qopenmp -qoverride-limits -xhost -DNDEBUG -std=c++11 -W
#
# ARPACK_INC = -larpack
# ARPACK_LIB = -larpack
#
# LAPACK_INC = -DUSE_MKL -DMARRAY_BLAS_PROTOTYPED -DMKL_Complex16="std::complex<double>" -I${MKLROOT}/include
# LAPACK_LIB = -L${MKLROOT}/lib -Wl,-rpath,${MKLROOT}/lib -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl
#
# HDF5_INC = -DUSE_EZH5 -I/usr/local/opt/hdf5/include
# HDF5_LIB = -L/usr/local/opt/hdf5/lib -lhdf5
#
# TBLIS_INC = -DUSE_TBLIS -I/usr/local/include -I/usr/local/include/tblis
# TBLIS_LIB = -L/usr/local/lib -ltblis

# CC = g++-7
# CC_FLAGS = -O2 -g -march=native -fopenmp -ffast-math -std=c++11 -W
#
# ARPACK_INC = -larpack
# ARPACK_LIB = -larpack
#
# LAPACK_INC = -framework Accelerate
# LAPACK_LIB = -framework Accelerate
#
# HDF5_INC = -DUSE_EZH5 -I/usr/local/opt/hdf5/include
# HDF5_LIB = -L/usr/local/opt/hdf5/lib -lhdf5
#
# TBLIS_INC = -DUSE_TBLIS -I/usr/local/include -I/usr/local/include/tblis
# TBLIS_LIB = -L/usr/local/lib -ltblis


CC = g++
CC_FLAGS = -O3 -DNDEBUG -march=native -fopenmp -ffast-math -std=c++11 -W

ARPACK_INC = -larpack
ARPACK_LIB = -larpack

LAPACK_INC = -I/usr/include/openblas
LAPACK_LIB = -lpthread -L/usr/lib/openblas-base -lopenblas

HDF5_INC = -DUSE_EZH5 -I/usr/lib/x86_64-linux-gnu/hdf5/serial/include
HDF5_LIB = -L/usr/lib/x86_64-linux-gnu/hdf5/serial/lib -lhdf5

TBLIS_INC = -I/usr/local/include -I/usr/local/include/tblis
TBLIS_LIB = -L/usr/local/lib -ltblis

# File names
EXEC = run
SOURCES = $(wildcard ./main.cpp)
OBJECTS = $(SOURCES:.cpp=.o)

# Main target
$(EXEC): $(OBJECTS)
	$(CC) $(OBJECTS) $(CC_FLAGS) $(TBLIS_LIB) $(ARPACK_LIB) $(LAPACK_LIB) $(HDF5_LIB) -o $(EXEC)

# To obtain object files
%.o: %.cpp
	$(CC) -c $(CC_FLAGS) $(TBLIS_INC) $(ARPACK_INC) $(LAPACK_INC) $(HDF5_INC) $< -o $@

# To remove generated files
clean:
	rm -f $(EXEC) $(OBJECTS)
