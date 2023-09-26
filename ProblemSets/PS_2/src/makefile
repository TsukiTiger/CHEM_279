# /src/makefile

# Compiler and Compiler Flags
CXX = g++
CXXFLAGS = -std=c++11 -Wall

# Source and Object Files
SRC = PS_2.cpp
OBJ = ../lib/PS_2.o

all: $(OBJ)

$(OBJ): $(SRC)
	$(CXX) $(CXXFLAGS) -c $(SRC) -o $(OBJ)

clean:
	rm -f $(OBJ)
