# /test/makefile

# Compiler and Compiler Flags
CXX = g++
CXXFLAGS = -std=c++11 -Wall

# Executables
UNIT_TEST = ../bin/test_unit
INTEGRATED_TEST = ../bin/test_integrated

all: $(UNIT_TEST) $(INTEGRATED_TEST)

$(UNIT_TEST): test_unit.cpp ../src/PS_1.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

$(INTEGRATED_TEST): test_integrated.cpp ../src/PS_1.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

clean:
	rm -f $(UNIT_TEST) $(INTEGRATED_TEST)
