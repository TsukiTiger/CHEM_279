#include "molecule.h"
#include <iostream>
#include <cmath>

BasisFunction::BasisFunction(const std::string& AO_type,
                             const arma::rowvec& center,
                             const arma::vec& lmn,
                             const arma::vec& exponents,
                             const arma::vec& contraction_coeffs)
        : AO_type(AO_type),
          center(center),
          lmn(lmn),
          exponents(exponents),
          contraction_coeffs(contraction_coeffs)
{
    calcNormConstants();
}

void BasisFunction::calcNormConstants() {
    // Calculate the overlap integrals for each exponent with itself
    double selfOverlap_1 = overlapIntegral3D(center, center, exponents(0), exponents(0), lmn, lmn);
    double selfOverlap_2 = overlapIntegral3D(center, center, exponents(1), exponents(1), lmn, lmn);
    double selfOverlap_3 = overlapIntegral3D(center, center, exponents(2), exponents(2), lmn, lmn);

    // Calculate the normalization constants
    norm_constants = {1.0 / sqrt(selfOverlap_1), 1.0 / sqrt(selfOverlap_2), 1.0 / sqrt(selfOverlap_3)};

}

Molecule::Molecule(const std::string& filename) {
    // Open the file
    ifstream infile(filename);
    assert(infile.good());

    // Read in the number of atoms
    infile >> nAtoms;

    atomicNumbers.resize(nAtoms);
    coordinates.resize(nAtoms, 3);

    std::string line;
    getline(infile, line); // Skip the rest of the line

    for (int i = 0; i < nAtoms; i++) {
        // Get line and split it
        getline(infile, line);
        istringstream iss(line);
        iss >> atomicNumbers(i);
        // Loop over the three coordinates x, y, z
        for (int j = 0; j < 3; j++) {
            iss >> coordinates(i, j);
        }
    }

    // Close the file
    infile.close();

    // Count carbon and hydrogen atoms
    carbon_count = 0;
    hydrogen_count = 0;
    for (int i = 0; i < nAtoms; i++) {
        if (atomicNumbers(i) == 1) {
            hydrogen_count++;
        } else if (atomicNumbers(i) == 6) {
            carbon_count++;
        }
    }

    // Calculate the number of electrons and basis functions
    nElectrons = countElectrons();
    nBasisFunctions = countBasisFunctions();

    // Information on H and C basis functions
    H_exponents = {3.42525091, 0.62391373, 0.16885540};
    H_coefficients = {0.15432897, 0.53532814, 0.44463454};
    C_exponents = {2.94124940, 0.68348310, 0.22228990};
    C_2s_coefficients = {-0.09996723, 0.39951283, 0.70011547};
    C_2p_coefficients = {0.15591627, 0.60768372, 0.39195739};

    // Build the list of basic functions
    basisFunctionsList = buildBasisFunctionsList();
}

void Molecule::printMoleculeInfo() const {
    cout << "Number of atoms: " << nAtoms << endl;
    cout << "Atomic numbers: " << atomicNumbers.t() << endl;
    cout << "Coordinates:\n" << coordinates << endl;
    cout << "Carbon count: " << carbon_count << endl;
    cout << "Hydrogen count: " << hydrogen_count << endl;
    cout << "Number of electrons: " << nElectrons << endl;
    cout << "Number of basis functions: " << nBasisFunctions << endl;
}

int Molecule::countBasisFunctions() {
    return 4 * carbon_count + hydrogen_count;
}

int Molecule::countElectrons() {
    double electron_pairs = 2 * carbon_count + (hydrogen_count / 2);

    // Check if electron_pairs is not an integer
    if (electron_pairs != int(electron_pairs)) {
        cout << "Error: Number of electron pairs is not an integer." << endl;
        exit(1);
    } else {
        return int(electron_pairs);
    }
}

std::vector<BasisFunction> Molecule::buildBasisFunctionsList() {
    // Initialize the list of basis functions
    std::vector<BasisFunction> basisFunctionsList;

    // Loop over all atoms
    for (int i = 0; i < nAtoms; i++) {

        // If the atom is hydrogen, add a 1s basis function
        if (atomicNumbers(i) == 1) {
            class BasisFunction AO_1s = {"H-1s", coordinates.row(i), {0, 0, 0}, H_exponents, H_coefficients};
            basisFunctionsList.push_back(AO_1s);

            // If the atom is carbon, add 2s and 2p basis functions (4 in total)
        } else if (atomicNumbers(i) == 6) {
            class BasisFunction AO_2s = {"C-2s", coordinates.row(i), {0, 0, 0}, C_exponents, C_2s_coefficients};
            class BasisFunction AO_2px = {"C-2px", coordinates.row(i), {1, 0, 0}, C_exponents, C_2p_coefficients};
            class BasisFunction AO_2py = {"C-2py", coordinates.row(i), {0, 1, 0}, C_exponents, C_2p_coefficients};
            class BasisFunction AO_2pz = {"C-2pz", coordinates.row(i), {0, 0, 1}, C_exponents, C_2p_coefficients};
            basisFunctionsList.push_back(AO_2s);
            basisFunctionsList.push_back(AO_2px);
            basisFunctionsList.push_back(AO_2py);
            basisFunctionsList.push_back(AO_2pz);
        }
    }
    return basisFunctionsList;
}

double overlapIntegral1D(double alpha, double beta, double center_a, double center_b, int lA, int lB) {
    // Calculate the exponential prefactor and the associated square root
    double prefactor = exp(-alpha * beta * pow(center_a - center_b, 2) / (alpha + beta));
    prefactor *= sqrt(M_PI / (alpha + beta));

    // Calculate center_product
    double center_product = (alpha * center_a + beta * center_b) / (alpha + beta);

    // Double summation over the angular momentum combinations
    double sum = 0.0;
    for (int i = 0; i <= lA; i++) {
        for (int j = 0; j <= lB; j++) {
            // Only (i + j) even terms contribute
            if ((i + j) % 2 == 0) {
                sum += binomialCoef(lA, i) * binomialCoef(lB, j) * (doubleFactorial(i + j - 1)
                                                                    * pow(center_product - center_a, lA - i) * pow(center_product - center_b, lB - j))
                       / pow(2 * (alpha + beta), double(i + j) / 2);
            }
        }
    }
    double integral = prefactor * sum;
    return integral;
}

double overlapIntegral3D(const arma::rowvec& centerA,
                         const arma::rowvec& centerB,
                         double expA, double expB,
                         const arma::vec& lmnA,
                         const arma::vec& lmnB)
{
    // Calculate the overlap integral for each dimension
    double integral = overlapIntegral1D(expA, expB, centerA(0), centerB(0), lmnA(0), lmnB(0)) *
                      overlapIntegral1D(expA, expB, centerA(1), centerB(1), lmnA(1), lmnB(1)) *
                      overlapIntegral1D(expA, expB, centerA(2), centerB(2), lmnA(2), lmnB(2));

    return integral;
}

double calcContractedOverlap(const BasisFunction& basisFunction1, const BasisFunction& basisFunction2) {
    double contracted_overlap = 0.0;

    // Loop over all exponents
    for (int k = 0; k < 3; k++) {
        for (int l = 0; l < 3; l++) {
            double unnorm_overlap = overlapIntegral3D(basisFunction1.center, basisFunction2.center,
                                                      basisFunction1.exponents(k), basisFunction2.exponents(l),
                                                      basisFunction1.lmn, basisFunction2.lmn);

            contracted_overlap += basisFunction1.contraction_coeffs[k] * basisFunction2.contraction_coeffs[l] *
                                  basisFunction1.norm_constants[k] * basisFunction2.norm_constants[l] * unnorm_overlap;
        }
    }
    return contracted_overlap;
}

arma::mat calcOverlapMatrix(const Molecule& molecule) {
    // Initialize the overlap matrix
    mat overlap_matrix = arma::zeros<mat>(molecule.nBasisFunctions, molecule.nBasisFunctions);

    // Loop over all basis function combinations
    for (int i = 0; i < molecule.nBasisFunctions; i++) {
        for (int j = 0; j < molecule.nBasisFunctions; j++) {
            overlap_matrix(i, j) = calcContractedOverlap(molecule.basisFunctionsList[i], molecule.basisFunctionsList[j]);
        }
    }

    return overlap_matrix;
}
