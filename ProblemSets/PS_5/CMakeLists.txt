cmake_minimum_required(VERSION 3.26)
project(EHT_Final)

set(CMAKE_CXX_STANDARD 17)

add_executable(EHT_Final main.cpp
        src/BasisFunction.cpp
        inc/BasisFunction.h
        src/Molecule.cpp
        inc/Molecule.h
        src/Atom.cpp
        inc/Atom.h
        src/CNDO.cpp
        inc/CNDO.h
        src/Utils.cpp
        inc/Utils.h
        inc/BasisSetInfo.h
        test/test.cpp)
