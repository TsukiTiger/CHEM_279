# /src/makefile

# Compiler and Compiler Flags
CXX = g++
CXXFLAGS = -std=c++11 -Wall

# Source and Object Files
SRC = PS_1.cpp
OBJ = ../lib/PS_1.o

all: $(OBJ)

$(OBJ): $(SRC)
	$(CXX) $(CXXFLAGS) -c $(SRC) -o $(OBJ)

clean:
	rm -f $(OBJ)
