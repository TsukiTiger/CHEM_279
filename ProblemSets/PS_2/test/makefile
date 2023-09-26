# /test/makefile

# Compiler and Compiler Flags
CXX = g++
CXXFLAGS = -std=c++11 -Wall

# Executables
UNIT_TEST = ../bin/test

all: $(UNIT_TEST)

$(UNIT_TEST): test.cpp ../src/PS_2.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

clean:
	rm -f $(UNIT_TEST)
