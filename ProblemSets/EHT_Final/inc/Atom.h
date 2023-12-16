//
// Created by Chongye Feng on 12/6/23.
//
#ifndef ATOM_H
#define ATOM_H

#include "BasisFunction.h"
#include <vector>
#include <string>

// Forward declaration
class BasisFunction;

class Atom {
public:
    Atom(const std::string& name, int valenceElectrons);

    void addBasisFunction(const BasisFunction& basisFunction);

    // Member variables
    std::string name;
    int valenceElectrons;
    std::vector<BasisFunction> mAOs;
};

#endif // ATOM_H

