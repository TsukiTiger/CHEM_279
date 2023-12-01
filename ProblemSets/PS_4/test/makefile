# Compiler and Compiler Flags
CXX = g++
CXXFLAGS = -std=c++11 -Wall -I../inc -I../utils -I/Users/chongyefeng/Documents/MSSE/CHEM_279/armadillo-12.6.2/include
LDFLAGS = -lblas -llapack

# Object Files
OBJ = ../lib/molecule.o ../lib/CNDO.o

# Executables
TEST = ../bin/test

all: $(TEST)

$(TEST): test.cpp $(OBJ)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

clean:
	rm -f $(TEST)
