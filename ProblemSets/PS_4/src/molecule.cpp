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

double I2e_pG(const arma::rowvec &Ra, const arma::rowvec &Rb, double sigmaA, double sigmaB)
{
    double U = pow(M_PI*sigmaA, 1.5) * pow(M_PI*sigmaB, 1.5); //eq 3.8 and 3.11
    double V2= 1.0/(sigmaA+sigmaB); //eq 3.9
    double Rd = arma::norm(Ra-Rb, 2);
    if (Rd==0.0) {
        //use eq 3.15
        return U * 2 * sqrt(V2/M_PI);
    }
    //eq 3.14 need sqrt T
    double srT = sqrt(V2) * Rd;
    //eq 3.14
    double result = U / Rd * erf(srT);
    return result;
}

double eval_2eI_s(const BasisFunction &ao1, const BasisFunction &ao2) {
    // Check if both ao1 and ao2 are s orbitals (lmn vectors should be {0,0,0})

    if (!(arma::accu(ao1.lmn) == 0) or !(arma::accu(ao2.lmn) == 0)) {
        std::cerr << "Error: Only s orbitals are allowed for eval_2eI_s function." << std::endl;
        return 0.0;
    }

    // Get the contraction coefficients multiplied by the normalization constants
    // Assuming da and db are arma::vec
    arma::vec da(ao1.contraction_coeffs.size());
    arma::vec db(ao2.contraction_coeffs.size());

    for (size_t i = 0; i < da.size(); ++i) {
        da[i] = ao1.contraction_coeffs[i] * ao1.norm_constants[i];
        db[i] = ao2.contraction_coeffs[i] * ao2.norm_constants[i];
    }


    // Extract exponents and centers
    const arma::vec& alphaa = ao1.exponents;
    const arma::vec& alphab = ao2.exponents;
    const arma::rowvec& Ra = ao1.center;
    const arma::rowvec& Rb = ao2.center;

    // Assuming the length of alpha vectors is 3 for STO-3G basis set
    const int len = 3;

    double sum = 0.0;
    // Nested loops to compute the sum over all primitives in the contracted Gaussians
    for (int k1 = 0; k1 < len; ++k1) {
        for (int k2 = 0; k2 < len; ++k2) {
            double sigmaA = 1.0 / (alphaa(k1) + alphaa(k2)); // Pre-factor for AO1

            for (int j1 = 0; j1 < len; ++j1) {
                for (int j2 = 0; j2 < len; ++j2) {
                    double sigmaB = 1.0 / (alphab(j1) + alphab(j2)); // Pre-factor for AO2
                    double I2e = I2e_pG(Ra, Rb, sigmaA, sigmaB); // Two-electron integral calculation

                    // Accumulate the result with contraction and normalization factors
                    sum += da(k1) * da(k2) * db(j1) * db(j2) * I2e;
                }
            }
        }
    }

    return sum; // Return the computed two-electron integral (gamma value)
}


void Molecule::eval_gammamat(arma::mat &gamma_mat) {
    gamma_mat.zeros(nAtoms, nAtoms); // Initialize gamma matrix with zeros

    for (size_t i = 0; i < nAtoms; i++) {
        for (size_t j = 0; j < nAtoms; j++) {
            // Find the s orbital for atom i
            const BasisFunction *s_orbital_i = findSOrbital(atoms[i]);

            // Find the s orbital for atom j
            const BasisFunction *s_orbital_j = findSOrbital(atoms[j]);

            if (s_orbital_i && s_orbital_j) {
                // Calculate the two-electron integral
                gamma_mat(i, j) = eval_2eI_s(*s_orbital_i, *s_orbital_j);
            }
        }
    }
}

// Helper function to find the s orbital of an atom
const BasisFunction* Molecule::findSOrbital(const Atom &atom) {
    for (auto &orbital : atom.mAOs) {
        if (arma::accu(orbital.lmn) == 0) { // Check if lmn is {0, 0, 0}
            return &orbital;
        }
    }
    return nullptr; // No s orbital found
}

void BasisFunction::calcNormConstants() {
    // Calculate the overlap integrals for each exponent with itself

    double selfOverlap_1 = Molecule::overlapIntegral3D(center, center, exponents(0), exponents(0), lmn, lmn);
    double selfOverlap_2 = Molecule::overlapIntegral3D(center, center, exponents(1), exponents(1), lmn, lmn);
    double selfOverlap_3 = Molecule::overlapIntegral3D(center, center, exponents(2), exponents(2), lmn, lmn);

    // Calculate the normalization constants
    norm_constants = {1.0 / sqrt(selfOverlap_1), 1.0 / sqrt(selfOverlap_2), 1.0 / sqrt(selfOverlap_3)};
}

Molecule::Molecule(const std::string& filename) {
    // Build the info maps for the molecule object
    buildBasisMap();

    // Open the file
    ifstream infile(filename);
    assert(infile.good());

    // Define a map from element symbols to atomic numbers for the first 10 elements
    std::map<std::string, int> elementToNumber = {
            {"H", 1}, {"He", 2}, {"Li", 3}, {"Be", 4}, {"B", 5},
            {"C", 6}, {"N", 7}, {"O", 8}, {"F", 9}, {"Ne", 10}
    };

    // Read in the number of atoms
    infile >> nAtoms;
    infile >> charge;

    elements.resize(nAtoms);
    atomicNumbers.resize(nAtoms);
    coordinates.resize(nAtoms, 3);

    std::string line;
    getline(infile, line); // Skip the rest of the line

    for (int i = 0; i < nAtoms; i++) {
        // Get line and split it
        getline(infile, line);
        istringstream iss(line);
        iss >> elements[i]; // Read element symbol

        // Map element symbol to atomic number
        atomicNumbers(i) = elementToNumber[elements[i]];

        // Loop over the three coordinates x, y, z
        for (int j = 0; j < 3; j++) {
            iss >> coordinates(i, j);
        }
    }

    // Close the file
    infile.close();

    // Calculate the number of electrons and basis functions
    countElectrons();
    countBasisFunctions();

    // Build the list of basic functions
    buildAtomsAndBasisFunctions();
}

void Molecule::printMoleculeInfo() const {
    cout << "Number of atoms: " << nAtoms << endl;
    cout << "Atomic numbers: " << atomicNumbers.t() << endl;
    cout << "Coordinates:\n" << coordinates << endl;
    cout << "Carbon count: " << carbon_count << endl;
    cout << "Hydrogen count: " << hydrogen_count << endl;
    cout << "Number of electron pairs: " << nElectrons << endl;
    cout << "Number of basis functions: " << nBasisFunctions << endl;
}

void Molecule::countBasisFunctions() {
    nBasisFunctions = 0; // Reset the count

    for (int i = 0; i < nAtoms; i++) {
        if (atomicNumbers(i) == 1) {  // Hydrogen
            nBasisFunctions += 1;
        } else {
            // For other elements like Carbon, Nitrogen, Oxygen, and Fluorine
            nBasisFunctions += 4;
        }
    }
}

// Which is valence electron
void Molecule::countElectrons() {
    nElectrons = 0; // Reset the count

    for (int i = 0; i < nAtoms; i++) {
        int atomicNumber = atomicNumbers(i);
        if (atomicNumber == 1) {  // Hydrogen
            nElectrons += 1;  // Hydrogen has 1 valence electron
        } else {
            // For elements like Carbon, Nitrogen, Oxygen, and Fluorine
            nElectrons += (atomicNumber - 2);  // Subtract 2 to get valence electrons
        }
    }

    // Adjust for molecular charge
    nElectrons -= charge;
}

void Molecule::buildAtomsAndBasisFunctions() {
    atoms.clear();
    for (int i = 0; i < nAtoms; i++) {
        std::string element = elements[i];
        Atom atom(element, atomicNumbers[i]); // Using atomicNumbers for valence electrons

        // Fetch the basis set information
        const BasisSetInfo& info = basisSetMap.at(element);

        // Add 's' orbitals
        std::string AO_label = element + "-2s";
        arma::rowvec center = coordinates.row(i);
        arma::vec lmn = {0, 0, 0}; // 's' orbital
        arma::vec exponents = info.exponents;
        arma::vec contraction_coeffs = info.s_coefficients;

        BasisFunction s_orbital(AO_label, center, lmn, exponents, contraction_coeffs);

        atom.addBasisFunction(s_orbital);

        basisFunctionsList.push_back(s_orbital);

        // Add 'p' orbitals if available
        if (!info.p_coefficients.empty()) {
            std::string AO_label = element + "-2p";
            arma::rowvec center = coordinates.row(i);

            // 'p' orbitals in x, y, and z directions
            arma::vec lmn_px = {1, 0, 0}; // 2px
            arma::vec lmn_py = {0, 1, 0}; // 2py
            arma::vec lmn_pz = {0, 0, 1}; // 2pz

            arma::vec exponents_p = {info.exponents};
            arma::vec contraction_coeffs_p = {info.p_coefficients};

            BasisFunction p_orbital_x(AO_label, center, lmn_px, exponents_p, contraction_coeffs_p);
            BasisFunction p_orbital_y(AO_label, center, lmn_py, exponents_p, contraction_coeffs_p);
            BasisFunction p_orbital_z(AO_label, center, lmn_pz, exponents_p, contraction_coeffs_p);

            atom.addBasisFunction(p_orbital_x);
            atom.addBasisFunction(p_orbital_y);
            atom.addBasisFunction(p_orbital_z);

            basisFunctionsList.push_back(p_orbital_x);
            basisFunctionsList.push_back(p_orbital_y);
            basisFunctionsList.push_back(p_orbital_z);
        }
        atoms.push_back(atom);
    }
}

double Molecule::overlapIntegral1D(double alpha, double beta, double center_a, double center_b, int lA, int lB) {
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

double Molecule::overlapIntegral3D(const arma::rowvec& centerA,
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

double Molecule::calcContractedOverlap(const BasisFunction& basisFunction1, const BasisFunction& basisFunction2) {
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

arma::mat Molecule::calcOverlapMatrix() {
    // Initialize the overlap matrix
    arma::mat overlap_matrix = arma::zeros<arma::mat>(nBasisFunctions, nBasisFunctions);

    std::cout << "Number of basis functions in list: " << basisFunctionsList.size() << std::endl;

    // Loop over all basis function combinations
    for (int i = 0; i < nBasisFunctions; i++) {
        for (int j = 0; j < nBasisFunctions; j++) {
            overlap_matrix(i, j) = calcContractedOverlap(basisFunctionsList[i], basisFunctionsList[j]);
        }
    }

    return overlap_matrix;
}

void Molecule::buildBasisMap() {
    basisSetMap["H"] = {
            {3.42525091, 0.62391373, 0.16885540},
            {0.15432897, 0.53532814, 0.44463454},
            {} // Hydrogen has no 'p' orbitals in STO-3G
    };

    basisSetMap["C"] = {
            {2.94124940, 0.68348310, 0.22228990},
            {-0.09996723, 0.39951283, 0.70011547},
            {0.15591627, 0.60768372, 0.39195739}
    };

    basisSetMap["N"] = {
            {3.78045590, 0.87849660, 0.28571440},
            {-0.09996723, 0.39951283, 0.70011547},
            {0.15591627, 0.60768372, 0.39195739}
    };

    basisSetMap["O"] = {
            {5.03315130, 1.16959610, 0.38038900},
            {-0.09996723, 0.39951283, 0.70011547},
            {0.15591627, 0.60768372, 0.39195739}
    };

    basisSetMap["F"] = {
            {6.46480320, 1.50228120, 0.48858850},
            {-0.09996723, 0.39951283, 0.70011547},
            {0.15591627, 0.60768372, 0.39195739}
    };
}

