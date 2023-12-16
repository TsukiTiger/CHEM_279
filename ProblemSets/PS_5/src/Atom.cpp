//
// Created by Chongye Feng on 12/6/23.
//

#include "../inc/Atom.h"


Atom::Atom(const std::string& name, int valenceElectrons)
        : name(name), valenceElectrons(valenceElectrons) {}

void Atom::addBasisFunction(const BasisFunction& basisFunction) {
    mAOs.push_back(basisFunction);
}
