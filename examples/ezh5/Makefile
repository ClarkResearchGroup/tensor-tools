# CC = g++
# CC_FLAGS = -I /usr/lib/x86_64-linux-gnu/hdf5/serial/include -L /usr/lib/x86_64-linux-gnu/hdf5/serial/lib -lhdf5 -std=c++11 -W

CC = icpc
CC_FLAGS = -O3 -qopenmp -xhost -DUSE_MKL -mkl=parallel -DUSE_EZH5 -I/usr/local/opt/hdf5/include -L/usr/local/opt/hdf5/lib -lhdf5 -std=c++11 -W

# File names
EXEC = run
SOURCES = $(wildcard ./main.cpp)
OBJECTS = $(SOURCES:.cpp=.o)

# Main target
$(EXEC): $(OBJECTS)
	$(CC) $(OBJECTS) $(CC_FLAGS) -o $(EXEC)

# To obtain object files
%.o: %.cpp
	$(CC) -c $(CC_FLAGS) $< -o $@

# To remove generated files
clean:
	rm -f $(EXEC) $(OBJECTS)
