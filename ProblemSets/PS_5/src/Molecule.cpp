//
// Created by Chongye Feng on 12/6/23.
//

#include "../inc/Molecule.h"

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
                sum += binomialCoef(lA, i)
                        * binomialCoef(lB, j)
                        * (doubleFactorial(i + j - 1)
                        * pow(center_product - center_a, lA - i)
                        * pow(center_product - center_b, lB - j))
                       / pow(2 * (alpha + beta), double(i + j) / 2);
            }
        }
    }
    double integral = prefactor * sum;
    return integral;
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

double Molecule::calcContractedOverlap(const BasisFunction& basisFunction1, const BasisFunction& basisFunction2) {
    double contracted_overlap = 0.0;

    // Loop over all exponents
    for (int k = 0; k < 3; k++) {
        for (int l = 0; l < 3; l++) {
            double unnorm_overlap = overlapIntegral3D(basisFunction1.getCenter(), basisFunction2.getCenter(),
                                                      basisFunction1.getExponents()(k), basisFunction2.getExponents()(l),
                                                      basisFunction1.getLmn(), basisFunction2.getLmn());

            contracted_overlap += basisFunction1.getContractionCoeffs()[k] * basisFunction2.getContractionCoeffs()[l] *
                                  basisFunction1.getNormConstants()[k] * basisFunction2.getNormConstants()[l] * unnorm_overlap;
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

const BasisFunction* Molecule::findSOrbital(const Atom &atom) {
    for (auto &orbital : atom.mAOs) {
        if (arma::accu(orbital.getLmn()) == 0) { // Check if lmn is {0, 0, 0}
            return &orbital;
        }
    }
    return nullptr; // No s orbital found
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
                gamma_mat(i, j) = eval_2eI_s( *s_orbital_i,  *s_orbital_j);
            }
        }
    }
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
    int ith_ao = 0;
    for (int i = 0; i < nAtoms; i++) {
        std::string element = elements[i];
        Atom atom(element, atomicNumbers[i]); // Using atomicNumbers for valence electrons

        // Fetch the basis set information
        const BasisSetInfo& info = basisSetMap.at(element);

        std::string AO_label;

        // Add 's' orbitals
        AO_label = element + "-2s";
        arma::rowvec center = coordinates.row(i);
        arma::vec lmn = {0, 0, 0}; // 's' orbital
        arma::vec exponents_s = info.exponents;
        arma::vec contraction_coeffs_s = info.s_coefficients;

        BasisFunction s_orbital(AO_label, center, lmn, exponents_s, contraction_coeffs_s);

        atom.addBasisFunction(s_orbital);

        basisFunctionsList.push_back(s_orbital);

        ao_map[ith_ao] = i;
        ith_ao += 1;

        // Add 'p' orbitals if available
        if (!info.p_coefficients.empty()) {
            AO_label = element + "-2p";
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
            ao_map[ith_ao] = i;
            ith_ao += 1;
            basisFunctionsList.push_back(p_orbital_y);
            ao_map[ith_ao] = i;
            ith_ao += 1;
            basisFunctionsList.push_back(p_orbital_z);
            ao_map[ith_ao] = i;
            ith_ao += 1;
        }
        atoms.push_back(atom);
    }
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
