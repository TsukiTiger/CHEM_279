# Compiler and Compiler Flags
CXX = g++
CXXFLAGS = -std=c++11 -Wall -O2 -I../inc/ -I/Users/chongyefeng/Documents/MSSE/CHEM_279/armadillo-12.6.2/include -DARMA_DONT_USE_WRAPPER
LDFLAGS = -framework Accelerate

# Source, Object and Header Files
SRC =  BasisFunction.cpp Utils.cpp Atom.cpp Molecule.cpp CNDO.cpp
OBJ = $(SRC:%.cpp=../lib/%.o)
HEADERS = $(wildcard ../inc/*.h)

all: $(OBJ)

../lib/%.o: %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f ../lib/*.o