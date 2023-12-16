//
// Created by Chongye Feng on 12/6/23.
//

#include "../inc/CNDO.h"

#include <iostream>
#include <vector>
#include <armadillo>
#include <unordered_map>
#include <string>

using namespace std;
using namespace arma;

void CNDO::initializeCNDOParams() {
    params["H"] = { {{"s", 7.176}}, 9.0 }; // H: Is + As and -β
    params["C"] = { {{"s", 14.051}, {"p", 5.572}}, 21.0 }; // C: (Is + As), (Ip + Ap) and -β
    params["N"] = { {{"s", 19.316}, {"p", 7.275}}, 25.0 }; // N: (Is + As), (Ip + Ap) and -β
    params["O"] = { {{"s", 25.390}, {"p", 9.111}}, 31.0 }; // O: (Is + As), (Ip + Ap) and -β
    params["F"] = { {{"s", 32.272}, {"p", 11.080}}, 39.0 }; // F: (Is + As), (Ip + Ap) and -β
}

void CNDO::rParamsInit() {
    // Initialize repulsion parameters for each atom pair (unit: eV)
    // Example values from Table 1 (actual values should be filled in)
    rparams["HH"] = {2.823,12.612, -0.0791,2.279};
    rparams["HC"] = {2.831,99.370, -0.0340,2.843};
    rparams["CC"] = {3.401,658.659, 0.0312,3.044};
    // Add more pairs as needed
}


CNDO::CNDO(Molecule &mol_i, int max_it, double tolerance)
        : mol(mol_i), max_iter(max_it), tol(tolerance) {
    // Initialize dimensions and matrices
    dim = mol.nBasisFunctions;

    Pa = arma::zeros<mat>(dim, dim);
    Pb = arma::zeros<mat>(dim, dim);
    Ca.set_size(dim, dim);
    Cb.set_size(dim, dim);
    Ea.set_size(dim);
    Eb.set_size(dim);
    H_core.set_size(dim, dim);

    // Populate the overlap matrix
    S = mol.calcOverlapMatrix();

    // Initialize and populate the gamma matrix
    int num_atoms = mol.nAtoms;
    gamma.set_size(num_atoms, num_atoms);

    // need to figure out what is the gamma? Is it gamma_AB?
    mol.eval_gammamat(gamma);
    gamma *= HARTREE_TO_EV; // Convert to electron volts


    // This Part is the
    // Calculate nuclear repulsion energy : Vnuc
    Ec = 0.0;
    for (int k = 0; k < static_cast<int>(mol.nAtoms); k++) {
        for (int j = k+1; j < static_cast<int>(mol.nAtoms); j++) {
            rowvec Ra = mol.atoms[k].mAOs[0].getCenter(); // Position of nucleus of atom k
            rowvec Rb = mol.atoms[j].mAOs[0].getCenter(); // Position of nucleus of atom j
            double Rd = arma::norm(Ra - Rb, 2);
            double ZA = mol.atomicNumbers[k];
            double ZB = mol.atomicNumbers[j];
            if(ZA > 2){
                ZA -= 2;
            }
            if(ZB > 2){
                ZB -= 2;
            }
            Ec += ZA * ZB / Rd; // Use atomic numbers for ZA and ZB
        }
    }
    Ec *= HARTREE_TO_EV; // Convert to electron volts

    // Determine the number of alpha and beta electrons
    q = mol.nElectrons / 2;
    p = mol.nElectrons - q;
    cout << "p = " << p << ", q = " << q << endl;
    // Given a way to allow the user to modify the p and q
    // Let the user input multiplicity

    // Initialize electronic and total energies
    Ee = 0.0;
    Etotal = 0.0;

    scfConverged = false;

    initializeCNDOParams();
    rParamsInit();

    initialize();
}


void CNDO::initialize() {
    // Initialize Fock matrices with appropriate dimensions
    F_alpha.set_size(dim, dim);
    F_beta.set_size(dim, dim);

    // Initialize Fock matrices to zero
    F_alpha.zeros();
    F_beta.zeros();

    // Initialize auxiliary matrices Ga and Gb with the same dimensions
    Ga.set_size(dim, dim);
    Gb.set_size(dim, dim);

    // Initialize Ga and Gb to zero
    Ga.zeros();
    Gb.zeros();

    // Build the core Hamiltonian matrix
    buildCoreHamiltonian(false);

    // Initialize the alpha and beta density matrices
    Pa.zeros(dim, dim);
    Pb.zeros(dim, dim);
}

void CNDO::buildCoreHamiltonian(bool repulsionCorrection) {
    H_core.zeros(dim, dim);
    size_t k_AO = 0;
    for (size_t k = 0; k < mol.nAtoms; k++) {
        const Atom &A_atom = mol.atoms[k];
        double ZA = static_cast<double>(mol.atomicNumbers(k));
        if (ZA > 2) {
            ZA -= 2; // counting only the valence
        }
        double gammaAA = gamma(k, k);

        for (const auto &ao_A : A_atom.mAOs) {
            string aoType(1, ao_A.getAoType().back());
            double IA = params[A_atom.name].ionizationAffinity[aoType];

            double sumZBgammaAB = 0.0;
            size_t j_AO = 0;
            for (size_t j = 0; j < mol.nAtoms; j++) {
                const Atom &B_atom = mol.atoms[j];
                double ZB = static_cast<double>(mol.atomicNumbers(j));
                if (ZB > 2) {
                    ZB -= 2; // counting only the valence
                }

                double beta_A = params[A_atom.name].beta;
                double beta_B = params[B_atom.name].beta;
                double beta_avg = -(beta_A + beta_B) / 2.0;

                if (k != j) {
                    sumZBgammaAB += ZB * gamma(k, j);
                }

                for (const auto &ao_B : B_atom.mAOs) {
                    if (k_AO != j_AO && k_AO < dim && j_AO < dim) {
                        H_core(k_AO, j_AO) = beta_avg * S(k_AO, j_AO);
                    }
                    j_AO++;
                }
            }

            if (k_AO < dim) {
                H_core(k_AO, k_AO) = -(IA) - (ZA - 0.5) * gammaAA - sumZBgammaAB;
            }
            k_AO++;
        }
    }
    cout << "Core Hamiltonian built successfully." << endl;
    if (repulsionCorrection){updateRepulsionHamiltonian();}
}

double CNDO::calculateRepulsionCorrection(double RAB, const std::string& atomPair) {
    // Retrieve parameters from the rparams map
    RepulsionParams params = rparams[atomPair];

    // Apply the repulsion correction formula
    double term1 = params.gamma * exp(-params.alpha * RAB);
    double term2 = params.omega * exp(-6 * pow((RAB - params.r), 2));
    return term1 + term2;
}

void CNDO::updateRepulsionHamiltonian(){
    // Eq (6) of the EHT paper
    // Correction of the repulsion
    for (size_t i = 0; i < mol.nAtoms; i++) {
        double repulsion_correction = 0.0;
        for (size_t j = 0; j < mol.nAtoms; j++) {
            if (i != j) {
                string atomPair = mol.atoms[i].name + mol.atoms[j].name;
                if (atomPair == "CH") {atomPair = "HC";}
                double distance = arma::norm(mol.atoms[i].mAOs[0].getCenter() - mol.atoms[j].mAOs[0].getCenter(), 2);
                repulsion_correction += calculateRepulsionCorrection(distance, atomPair);
            }
        }
        H_core(i, i) += repulsion_correction; // Update diagonal element
    }
}

void CNDO::calculateFockMatrix() {
    size_t u_AO = 0; // Index for atomic orbitals
    for (int A = 0; A < static_cast<int>(mol.nAtoms); A++) {
        Atom &atom_A = mol.atoms[A];
        string elementType = atom_A.name; // Element type of atom A
        double ZA = static_cast<double>(mol.atomicNumbers(A)); // Atomic number of atom A

        for (auto &ao_A : atom_A.mAOs) {
            string aoType = string(1, ao_A.getAoType().back()); // Get the last character of AO_type
            // Type of atomic orbital (s, p, etc.)
            double IA = params[elementType].ionizationAffinity[aoType]; // 1/2 * (Iu + Au)
            double beta_A = params[elementType].beta; // -beta

            double p_tot_AA_alpha = Pa(u_AO, u_AO); // Total alpha electron density at orbital μ
            double p_tot_AA_beta = Pb(u_AO, u_AO); // Total beta electron density at orbital μ
            double gamma_AA = gamma(A, A); // γ_AA value

            // Diagonal elements of the Fock matrix for alpha and beta electrons
            F_alpha(u_AO, u_AO) = -IA + ((p_tot_AA_alpha - ZA) - (Pa(u_AO, u_AO) - 0.5)) * gamma_AA;
            F_beta(u_AO, u_AO) = -IA + ((p_tot_AA_beta - ZA) - (Pb(u_AO, u_AO) - 0.5)) * gamma_AA;

            // Off-diagonal elements for B != A
            for (int B = 0; B < static_cast<int>(mol.nAtoms); B++) {
                if (B != A) {
                    double gamma_AB = gamma(A, B); // γ_AB value
                    for (int v_AO = 0; v_AO < dim; ++v_AO) {
                        if (u_AO != v_AO) {
                            double beta_B = params[mol.atoms[B].name].beta; // Beta parameter for atom B
                            double beta_avg = -(beta_A + beta_B) / 2.0; // Average beta for atoms A and B

                            F_alpha(u_AO, v_AO) = beta_avg * S(u_AO, v_AO) - Pa(u_AO, v_AO) * gamma_AB;
                            F_beta(u_AO, v_AO) = beta_avg * S(u_AO, v_AO) - Pb(u_AO, v_AO) * gamma_AB;
                        }
                    }
                }
            }
            u_AO++;
        }
    }
}


double CNDO::calculateTotalEnergy() {
    // Validation (optional, based on your overall program structure)
    if (Pa.n_rows != Pb.n_rows || Pa.n_cols != Pb.n_cols ||
        H_core.n_rows != Pa.n_rows || H_core.n_cols != Pa.n_cols ||
        F_alpha.n_rows != Pa.n_rows || F_alpha.n_cols != Pa.n_cols ||
        F_beta.n_rows != Pa.n_rows || F_beta.n_cols != Pa.n_cols) {
        throw std::runtime_error("Matrix dimension mismatch in calculateTotalEnergy");
    }

    // Ensure Ec is calculated and valid
    if (std::isnan(Ec) || std::isinf(Ec)) {
        throw std::runtime_error("Invalid nuclear repulsion energy (Ec)");
    }

    // Electronic energy calculation
    arma::mat Ptotal = Pa + Pb;
    Ee = arma::dot(Pa, Ga) / 2.0 + arma::dot(Pb, Gb) / 2.0 + arma::dot(Ptotal, H_core);
    Etotal = Ee + Ec;

    // Print nuclear repulsion energy
    std::cout << "Nuclear Repulsion Energy is " << Ec << " eV." << std::endl;
    std::cout << "Electron Energy is " << Ee << " eV." << std::endl;

    return Etotal;
}

void CNDO::runSCF() {
    // Step 1: Initial Guess
    Pa.zeros();
    Pb.zeros();

    // Initialize auxiliary G matrices
    Ga.set_size(dim, dim);
    Gb.set_size(dim, dim);

    // Init F matrix
    F_alpha = H_core + Ga;
    F_beta = H_core + Gb;

    // Iterate until convergence or max iterations
    for (size_t iter = 0; iter < max_iter; ++iter) {
        cout << "SCF Iteration: " << iter + 1 << endl;

        // Update G matrices based on current Pa and Pb
        updateG();

        // Step 2: Build Fock matrices
        F_alpha = H_core + Ga;
        F_beta = H_core + Gb;

        // Step 3: Solve eigenvalue problems
        arma::eig_sym(Ea, Ca, F_alpha);
        arma::eig_sym(Eb, Cb, F_beta);

        // Step 4: Copy old density matrices
        arma::mat Pa_old = Pa;
        arma::mat Pb_old = Pb;

        // Step 5: Assemble new density matrices
        Pa = Ca.cols(0, p - 1) * trans(Ca.cols(0, p - 1));
        if (q > 0) {
            Pb = Cb.cols(0, q - 1) * trans(Cb.cols(0, q - 1));
        } else {
            Pb.zeros();
        }

        // Check for convergence
        if (arma::approx_equal(Pa, Pa_old, "absdiff", tol) &&
            arma::approx_equal(Pb, Pb_old, "absdiff", tol)) {
            cout << "SCF Converged in " << iter + 1 << " iterations." << endl;
            scfConverged = true;
            break;
        }
    }

    // Calculate total energy after convergence
    Etotal = calculateTotalEnergy();
    cout << "Total Energy: " << Etotal << " eV" << endl;

    calculateGradient();
    printGradients();
}

void CNDO::updateG() {
    Ga.zeros();
    Gb.zeros();

    arma::vec P_t = arma::zeros(mol.nAtoms);
    size_t k_AO = 0;
    // Calculate total density for each atom
    for (size_t k = 0; k < mol.nAtoms; k++) {
        for (auto &ao : mol.atoms[k].mAOs) {
            P_t(k) += Pa(k_AO, k_AO) + Pb(k_AO, k_AO);
            k_AO++;
        }
    }

    k_AO = 0; // Resetting AO index for main loop
    for (size_t k = 0; k < mol.nAtoms; k++) {
        double gammaAA = gamma(k, k);
        for (auto &ao_A : mol.atoms[k].mAOs) {
            Ga(k_AO, k_AO) = (P_t(k) - Pa(k_AO, k_AO)) * gammaAA;
            Gb(k_AO, k_AO) = (P_t(k) - Pb(k_AO, k_AO)) * gammaAA;

            size_t j_AO = 0;
            for (size_t j = 0; j < mol.nAtoms; j++) {
                double gammaAB = gamma(k, j);
                if (k != j){
                    Ga(k_AO, k_AO) += P_t[j] * gammaAB;
                    Gb(k_AO, k_AO) += P_t[j] * gammaAB;
                }
                for (auto &ao_B : mol.atoms[j].mAOs) {
                    if (k_AO != j_AO) {
                        Ga(k_AO, j_AO) = -gammaAB * Pa(k_AO, j_AO);
                        Gb(k_AO, j_AO) = -gammaAB * Pb(k_AO, j_AO);
                    }
                    j_AO++;
                }
            }
            k_AO++; // Increment AO index
        }
    }
//    Ga.print("Ga: ");
//    Gb.print("Gb: ");
}


void CNDO::printResults() const {
    cout << "Gamma Matrix:" << endl;
    gamma.print();

    cout << "Overlap Matrix:" << endl;
    S.print();

    cout << "p = " << p << ", q = " << q << endl;

    cout << "Core Hamiltonian (H_core):" << endl;
    H_core.print();

    cout << "Fock Matrix (Alpha):" << endl;
    F_alpha.print();

    cout << "Fock Matrix (Beta):" << endl;
    F_beta.print();

    cout << "Alpha Orbital Coefficients (Ca):" << endl;
    Ca.print();

    cout << "Beta Orbital Coefficients (Cb):" << endl;
    Cb.print();

    cout << "Alpha Orbital Energies (Ea):" << endl;
    Ea.print();

    cout << "Beta Orbital Energies (Eb):" << endl;
    Eb.print();

    cout << "Alpha Density Matrix (Pa):" << endl;
    Pa.print();

    cout << "Beta Density Matrix (Pb):" << endl;
    Pb.print();

    cout << "Nuclear Repulsion Energy: " << Ec << " eV" << endl;
    cout << "Electronic Energy: " << Ee << " eV" << endl;
    cout << "Total Energy: " << Etotal << " eV" << endl;
}

void CNDO::printGradients() const {
    std::cout << "Nuclear Gradients: " << std::endl;
    n_gradient.print();
    std::cout << "Electronic Gradients: " << std::endl;
    e_gradient.print();
    std::cout << "Gradients: " << std::endl;
    gradient.print();
}

void CNDO::calculateGradient() {
    // Make sure the SCF has converged and you have the final density matrices
    if (!scfConverged) {
        std::cerr << "SCF not converged. Gradient cannot be computed." << std::endl;
        return;
    }

    // Calculate the derivative of the overlap integral
    calculateOV_RA();
    cout << "OV_RA Done" << endl;

    // Calculate the derivative of the gamma matrix
    calculateGamma_RA();
    cout << "gamma_RA Done" << endl;

    // Assemble the x matrix from dS_RA
    assembleXMatrix();
    cout << "x Done" << endl;

    // Assemble the y matrix from dGamma_RA
    assembleYMatrix();
    cout << "y Done" << endl;

    // Calculate the nuclear repulsion gradient
    n_gradient = calculateNuclearRepulsionGradient();
    cout << "n_g Done" << endl;

    // Calculate the electronic gradient
    e_gradient = calculateElectronicGradient();
    cout << "e_g Done" << endl;

    // Combine all the contributions to get the full gradient
    gradient.set_size(3, mol.nAtoms); // Ensure the gradient has the correct size
    gradient.zeros(); // Initialize to zero

    // Loop over each atom and sum the contributions to the gradient
    for (int atomIdx = 0; atomIdx < mol.nAtoms; ++atomIdx) {
        gradient.col(atomIdx) = e_gradient.col(atomIdx) + n_gradient.col(atomIdx);
    }
    cout << "g Done" << endl;
}

void CNDO::calculateOV_RA() {
    dS_RA.set_size(mol.nBasisFunctions, mol.nBasisFunctions);

    for (int m = 0; m < mol.nBasisFunctions; m++){
        for (int n = 0; n < mol.nBasisFunctions; n++){
            if (m != n){
                dS_RA(m, n) = OV_vec(mol.basisFunctionsList[m], mol.basisFunctionsList[n]);
            }
            else{
                dS_RA(m, n) = {0, 0, 0};
            }
        }
    }
}


vec CNDO::OV_vec(const BasisFunction& a, const BasisFunction& b)
{
    double dx = 0;
    double dy = 0;
    double dz = 0;

    vec result = {0, 0, 0};

    //x dim
    for (int k = 0; k < 3; k++){
        for (int l = 0; l < 3; l++){
            dx += a.getContractionCoeffs()[k] * b.getContractionCoeffs()[l]
                   * a.getNormConstants()[k] * b.getNormConstants()[l]
                   * d_s_AB(a.getExponents()[k], b.getExponents()[l], a.getCenter()[k], b.getCenter()[l], a.getLmn()[0], b.getLmn()[0])
                   * Molecule::overlapIntegral1D(a.getExponents()[k], b.getExponents()[l], a.getCenter()[k], b.getCenter()[l], a.getLmn()[1], b.getLmn()[1])
                   * Molecule::overlapIntegral1D(a.getExponents()[k], b.getExponents()[l], a.getCenter()[k], b.getCenter()[l], a.getLmn()[2], b.getLmn()[2]);
        }
    }

    // y dim
    for (int k = 0; k < 3; k++)
    {
        for (int l = 0; l < 3; l++)
        {
            dy += a.getContractionCoeffs()[k] * b.getContractionCoeffs()[l]
                  * a.getNormConstants()[k] * b.getNormConstants()[l]
                  * d_s_AB(a.getExponents()[k], b.getExponents()[l], a.getCenter()[k], b.getCenter()[l], a.getLmn()[1], b.getLmn()[1])
                  * Molecule::overlapIntegral1D(a.getExponents()[k], b.getExponents()[l], a.getCenter()[k], b.getCenter()[l], a.getLmn()[0], b.getLmn()[0])
                  * Molecule::overlapIntegral1D(a.getExponents()[k], b.getExponents()[l], a.getCenter()[k], b.getCenter()[l], a.getLmn()[2], b.getLmn()[2]);
        }
    }

    // z dim
    for (int k = 0; k < 3; k++)
    {
        for (int l = 0; l < 3; l++)
        {
            dz += a.getContractionCoeffs()[k] * b.getContractionCoeffs()[l]
                  * a.getNormConstants()[k] * b.getNormConstants()[l]
                  * d_s_AB(a.getExponents()[k], b.getExponents()[l], a.getCenter()[k], b.getCenter()[l], a.getLmn()[2], b.getLmn()[2])
                  * Molecule::overlapIntegral1D(a.getExponents()[k], b.getExponents()[l], a.getCenter()[k], b.getCenter()[l], a.getLmn()[1], b.getLmn()[1])
                  * Molecule::overlapIntegral1D(a.getExponents()[k], b.getExponents()[l], a.getCenter()[k], b.getCenter()[l], a.getLmn()[0], b.getLmn()[0]);
        }
    }

    result(0) = dx;
    result(1) = dy;
    result(2) = dz;

    return result;
}

double CNDO::d_s_AB(double exp_a, double exp_b, double center_a, double center_b, int lA, int lB){
    double dx = 0;
    double center = (exp_a * center_a + exp_b * center_b) / (exp_a + exp_b);

    if (lA == 0)
    {
        dx = 2 * exp_a * Molecule::overlapIntegral1D(exp_a, exp_b, center_a, center_b, lA+1, lB);
    }

    if (lA == 1)
    {
        dx = (-1.0) * Molecule::overlapIntegral1D(exp_a, exp_b, center_a, center_b, lA-1, lB)
                + 2 * exp_a * Molecule::overlapIntegral1D(exp_a, exp_b, center_a, center_b, lA+1, lB);
    }
    return dx;
}

vec CNDO::d00_vec(const BasisFunction &a, const BasisFunction &b, double sigma_a, double sigma_b) {

    double u_ab = std::pow(M_PI * sigma_a, 1.5) * std::pow(M_PI * sigma_b, 1.5);
    double rd = arma::norm(a.getCenter() - b.getCenter(), 2);
    double v2 = 1.0 / (sigma_a + sigma_b);
    double srT = std::sqrt(v2) * rd;
    double T = std::pow(srT, 2);

    arma::vec d_0_0(3);

    if (rd != 0.0){
        d_0_0 = (u_ab) * (1.0/std::pow(rd, 2)) * ((2*std::sqrt(v2) *exp(-T)) / std::sqrt(M_PI) - erf(srT)/rd) * (a.getCenter().t() - b.getCenter().t());
        return d_0_0;
    }
    else{
        d_0_0 = {0.0, 0.0, 0.0};
        return d_0_0;
    }
}

vec CNDO::Gamma_vec(const BasisFunction &a, const BasisFunction &b){
    vec d_a_coefficients(a.getContractionCoeffs().size());
    vec d_b_coefficients(b.getContractionCoeffs().size());

    vec alpha_a_coefficients = a.getExponents();
    vec alpha_b_coefficients = b.getExponents();


    //calculate d prime
    for (int i = 0; i < 3; i++)
    {
        d_a_coefficients(i) = a.getContractionCoeffs()[i] * a.getNormConstants()[i];
        d_b_coefficients(i) = b.getContractionCoeffs()[i] * b.getNormConstants()[i];
    }

    vec result(3);
    result.zeros();

    for (int k_1 = 0; k_1 < 3; k_1++)
    {
        for (int k_2 = 0; k_2 < 3; k_2++)
        {
            double sigma_a = 1.0 / (alpha_a_coefficients(k_1) + alpha_a_coefficients(k_2));

            for (int l_1 = 0; l_1 < 3; l_1++)
            {
                for (int l_2 = 0; l_2 < 3; l_2++)
                {
                    double sigma_b = 1.0 / (alpha_b_coefficients(l_1) + alpha_b_coefficients(l_2));
                    result += d_a_coefficients(k_1) * d_a_coefficients(k_2)
                              * d_b_coefficients(l_1) * d_b_coefficients(l_2)
                              * d00_vec(a, b, sigma_a, sigma_b);
                }
            }
        }
    }
    return result;
}

void CNDO::calculateGamma_RA() {
    dGamma_RA.set_size(mol.nAtoms, mol.nAtoms);

    for (int i = 0; i < mol.nAtoms; i++) {
        for (int j = 0; j < mol.nAtoms; j++) {
            const BasisFunction *a = mol.findSOrbital(mol.atoms[i]);
            const BasisFunction *b = mol.findSOrbital(mol.atoms[j]);
            if (a && b) {
                vec result = Gamma_vec(*a, *b);
                dGamma_RA(i, j) = result * HARTREE_TO_EV;
            }
        }
    }
}

void CNDO::assembleXMatrix() {
    xMatrix.set_size(mol.nBasisFunctions, mol.nBasisFunctions);

    for (int i_ao=0; i_ao<mol.nBasisFunctions; i_ao++){
        for (int j_ao=0; j_ao<mol.nBasisFunctions; j_ao++){
            string nameA = string(1, mol.basisFunctionsList[i_ao].getAoType().at(0));
            string nameB = string(1, mol.basisFunctionsList[j_ao].getAoType().at(0));
            double BA = params[nameA].beta;
            double BB = params[nameB].beta;

            xMatrix(i_ao, j_ao) = (BA + BB) * (Pa(i_ao, j_ao) + Pb(i_ao, j_ao));
        }
    }
}

void CNDO::assembleYMatrix() {
    yMatrix.set_size(mol.nAtoms, mol.nAtoms);


    for (int i_atom=0; i_atom<mol.nAtoms; i_atom++){
        for (int j_atom=0; j_atom<mol.nAtoms; j_atom++){
            double ZA = mol.atomicNumbers(i_atom);
            double ZB = mol.atomicNumbers(j_atom);

            if (ZA>2) {ZA-=2;}
            if (ZB>2) {ZB-=2;}

            double PA_tot = 0;
            double PB_tot = 0;
            double sum_part = 0;
            for (int i=0; i<mol.nBasisFunctions; i++){
                for (int j=0; j<mol.nBasisFunctions; j++){
                    if (mol.ao_map[i] == i_atom && mol.ao_map[j] == j_atom && i != j){
                        sum_part += Pa(i, j) * Pa(i, j) + Pb(i, j) * Pb(i, j);
                    }
                    if (mol.ao_map[i] == i_atom && i == j){
                        PA_tot += Pa(i, i) + Pb(i, i);
                    }
                    if (mol.ao_map[j] == j_atom && i == j){
                        PB_tot += Pa(j, j) + Pb(j, j);
                    }
                }
            }

            if(i_atom == j_atom) {sum_part = 0;}

            yMatrix(i_atom, j_atom) = PA_tot * PB_tot - ZA * PB_tot - ZB * PA_tot - sum_part;
        }
    }
}

mat CNDO::calculateNuclearRepulsionGradient() {
    n_gradient.set_size(3, mol.nAtoms);
    n_gradient.zeros();

    for (size_t i = 0; i < mol.nAtoms; ++i) {
        for (size_t j = 0; j < mol.nAtoms; ++j) {
            if (i != j) {
                rowvec pos_i = mol.atoms[i].mAOs[0].getCenter();
                rowvec pos_j = mol.atoms[j].mAOs[0].getCenter();
                rowvec r_vec = pos_i - pos_j; // Vector from j to i
                double r_norm = arma::norm(r_vec, 2);

                double ZA = (mol.atomicNumbers[i] > 2) ? mol.atomicNumbers[i] - 2 : mol.atomicNumbers[i];
                double ZB = (mol.atomicNumbers[j] > 2) ? mol.atomicNumbers[j] - 2 : mol.atomicNumbers[j];

                rowvec gradient_contribution = -ZA * ZB * (std::pow(r_norm, -3.0)) * r_vec * HARTREE_TO_EV;
                n_gradient.col(i) += gradient_contribution.t();
            }
        }
    }
    return n_gradient;
}

mat CNDO::calculateElectronicGradient(){
    e_gradient.set_size(3, mol.nAtoms);
    vec sum(3);
    for (int i_atom=0; i_atom<mol.nAtoms; i_atom++){
        sum = {0,0,0};
        for (int i=0; i<mol.nBasisFunctions; i++){
            for (int j=0; j<mol.nBasisFunctions; j++){
                if (mol.ao_map[i] == i_atom && mol.ao_map[j] != i_atom){
                    sum += xMatrix(i,j) * dS_RA(i,j);
                }
            }
        }
        for (int j_atom=0; j_atom<mol.nAtoms; j_atom++){
            if (i_atom != j_atom){
                sum += yMatrix(i_atom, j_atom) * dGamma_RA(i_atom, j_atom);
            }
        }
        e_gradient.col(i_atom) = sum;
    }
    return e_gradient;
}
