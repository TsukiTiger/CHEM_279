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

class Atom {
public:
    std::string name;                // Chemical name (e.g., "C" for Carbon)
    int valenceElectrons;            // Number of valence electrons
    std::vector<BasisFunction> mAOs; // Vector of associated basis functions (atomic orbitals)

    Atom(const std::string& name, int valenceElectrons)
            : name(name), valenceElectrons(valenceElectrons) {}

    void addBasisFunction(const BasisFunction& basisFunction) {
        mAOs.push_back(basisFunction);
    }
};

struct BasisSetInfo {
    std::vector<double> exponents;
    std::vector<double> s_coefficients; // For 's' orbitals
    std::vector<double> p_coefficients; // For 'p' orbitals
};

double I2e_pG(const arma::rowvec &Ra, const arma::rowvec &Rb, double sigmaA, double sigmaB);
double eval_2eI_s(BasisFunction &ao1, BasisFunction &ao2);

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

    int get_charge(){return charge;};
    const BasisFunction* findSOrbital(const Atom &atom);
    void eval_gammamat(arma::mat &gamma_mat);
private:
    void countBasisFunctions();
    void countElectrons();
    vector<BasisFunction> buildBasisFunctionsList();
    void buildAtomsAndBasisFunctions();
    int charge;
};

#endif // MOLECULE_H
