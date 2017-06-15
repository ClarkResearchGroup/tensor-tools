# CC = icpc
# CC_FLAGS = -O3 -qopenmp -xhost -DNDEBUG -std=c++11 -W
#
# LIB_ARPACK = -larpack
#
# LAPACK_INC = -DUSE_MKL -I${MKLROOT}/include
# LAPACK_LIB = -L${MKLROOT}/lib -Wl,-rpath,${MKLROOT}/lib -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl
#
# HDF5_INC = -DUSE_EZH5 -I/usr/local/opt/hdf5/include
# HDF5_LIB = -L/usr/local/opt/hdf5/lib -lhdf5
#
# TBLIS_INC = -I/usr/local/include -I/usr/local/include/tblis
# TBLIS_LIB = -L/usr/local/lib -ltblis


CC = g++
CC_FLAGS = -O3 -march=native -fopenmp -ffast-math -DNDEBUG -std=c++11 -W

LIB_ARPACK = -larpack

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
	$(CC) $(OBJECTS) $(CC_FLAGS) $(TBLIS_LIB) $(LAPACK_LIB) $(HDF5_LIB) -o $(EXEC)

# To obtain object files
%.o: %.cpp
	$(CC) -c $(CC_FLAGS) $(TBLIS_INC) $(LAPACK_INC) $(HDF5_INC) $< -o $@

# To remove generated files
clean:
	rm -f $(EXEC) $(OBJECTS)
