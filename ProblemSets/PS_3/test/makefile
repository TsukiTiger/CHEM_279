# Compiler and Compiler Flags
CXX = g++
CXXFLAGS = -std=c++11 -Wall -I../inc -I../utils -I/Users/chongyefeng/Documents/MSSE/CHEM_279/armadillo-12.6.2/include
LDFLAGS = -lblas -llapack

# Object Files
OBJ = ../lib/molecule.o ../lib/eht_matrices.o

# Executables
TEST = ../bin/test
ENERGY_DIFF = ../bin/energy_diff

all: $(TEST) $(ENERGY_DIFF)

$(TEST): test.cpp $(OBJ)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

$(ENERGY_DIFF): energy_diff.cpp $(OBJ)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

clean:
	rm -f $(TEST) $(ENERGY_DIFF)
