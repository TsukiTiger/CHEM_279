#ifndef MOLECULE_H
#define MOLECULE_H

#include <armadillo>
#include <vector>
#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <cassert>
#include "factorial.h"

using namespace arma;
using namespace std;

// Class for individual basis functions
class BasisFunction {
    friend class Molecule;
public:
    BasisFunction(const std::string& AO_type,
                  const arma::rowvec& center,
                  const arma::vec& lmn,
                  const arma::vec& exponents,
                  const arma::vec& contraction_coeffs);
    std::string AO_type;
    arma::rowvec center;
    arma::vec lmn;
    arma::vec exponents;
    arma::vec contraction_coeffs;
    std::vector<double> norm_constants;

private:
    void calcNormConstants();
};

double overlapIntegral3D(const arma::rowvec& centerA,
                         const arma::rowvec& centerB,
                         double expA, double expB,
                         const arma::vec& lmnA,
                         const arma::vec& lmnB);

// Class that contains all the information about a molecule
class Molecule {
    friend class BasisFunction;
public:
    Molecule(const std::string& filename);
    int nAtoms;
    vec atomicNumbers;
    mat coordinates;
    int carbon_count;
    int hydrogen_count;
    int nElectrons;
    int nBasisFunctions;
    vector<BasisFunction> basisFunctionsList;

    // Information on H and C basis functions
    vec H_exponents;
    vec H_coefficients;
    vec C_exponents;
    vec C_2s_coefficients;
    vec C_2p_coefficients;

    void printMoleculeInfo() const;

private:
    int countBasisFunctions();
    int countElectrons();
    vector<BasisFunction> buildBasisFunctionsList();
};

double overlapIntegral1D(int l1, int l2, double PA, double PB,
                         double alphaP, double normA, double normB);

double overlapIntegral3D(const BasisFunction& bf1, const BasisFunction& bf2,
                         const arma::mat& atomCoords);

double calcContractedOverlap(const BasisFunction& bf1, const BasisFunction& bf2,
                             const arma::mat& atomCoords);

arma::mat calcOverlapMatrix(const Molecule& mol);

#endif // MOLECULE_H
