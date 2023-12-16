//
// Created by Chongye Feng on 12/6/23.
//

#ifndef EHT_FINAL_MOLECULE_H
#define EHT_FINAL_MOLECULE_H

#include "Atom.h"
#include "BasisSetInfo.h"
#include "BasisFunction.h"
#include <armadillo>
#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include "Utils.h"


using namespace arma;
using namespace std;

class BasisFunction;
class Atom;

// Class that contains all the information about a molecule
class Molecule {
    friend class BasisFunction;
public:
    Molecule(const std::string& filename);
    int nAtoms;
    vec atomicNumbers;
    std::vector<std::string> elements;
    mat coordinates;
    int carbon_count;
    int hydrogen_count;
    int nElectrons;
    int nBasisFunctions;
    vector<BasisFunction> basisFunctionsList;
    vector<Atom> atoms;
    // Information on H and C basis functions
    vec H_exponents;
    vec H_coefficients;
    vec C_exponents;
    vec C_2s_coefficients;
    vec C_2p_coefficients;

    void printMoleculeInfo() const;

    static double overlapIntegral1D(double alpha, double beta, double center_a, double center_b, int lA, int lB);

    static double overlapIntegral3D(const arma::rowvec& centerA,
                                    const arma::rowvec& centerB,
                                    double expA, double expB,
                                    const arma::vec& lmnA,
                                    const arma::vec& lmnB);

    double calcContractedOverlap(const BasisFunction& basisFunction1, const BasisFunction& basisFunction2);

    arma::mat calcOverlapMatrix();

    std::map<std::string, BasisSetInfo> basisSetMap;
    void buildBasisMap();

    int get_charge() const{return charge;};
    const BasisFunction* findSOrbital(const Atom &atom);
    void eval_gammamat(arma::mat &gamma_mat);
    map<int, int> ao_map;
private:
    void countBasisFunctions();
    void countElectrons();
    void buildAtomsAndBasisFunctions();
    int charge;
};


#endif //EHT_FINAL_MOLECULE_H
