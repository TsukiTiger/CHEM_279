#include <iostream>
#include <vector>
#include <armadillo>
#include <unordered_map>
#include <string>

#include "CNDO.h"

using namespace std;
using namespace arma;

void CNDO::initializeCNDOParams() {
    params["H"] = { {{"s", 7.176}}, 9.0 }; // H: Is + As and -β
    params["C"] = { {{"s", 14.051}, {"p", 5.572}}, 21.0 }; // C: (Is + As), (Ip + Ap) and -β
    params["N"] = { {{"s", 19.316}, {"p", 7.275}}, 25.0 }; // N: (Is + As), (Ip + Ap) and -β
    params["O"] = { {{"s", 25.390}, {"p", 9.111}}, 31.0 }; // O: (Is + As), (Ip + Ap) and -β
    params["F"] = { {{"s", 32.272}, {"p", 11.080}}, 39.0 }; // F: (Is + As), (Ip + Ap) and -β
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

    // Calculate nuclear repulsion energy : Vnuc
    Ec = 0.0;
    for (size_t k = 0; k < mol.nAtoms; k++) {
        for (size_t j = k+1; j < mol.nAtoms; j++) {
            rowvec Ra = mol.atoms[k].mAOs[0].center; // Position of nucleus of atom k
            rowvec Rb = mol.atoms[j].mAOs[0].center; // Position of nucleus of atom j
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

    // Initialize electronic and total energies
    Ee = 0.0;
    Etotal = 0.0;

    initializeCNDOParams();
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
    buildCoreHamiltonian();

    // Initialize the alpha and beta density matrices
    Pa.zeros(dim, dim);
    Pb.zeros(dim, dim);
}

void CNDO::buildCoreHamiltonian() {
    H_core.zeros(dim, dim);
    size_t k_AO = 0;
    for (size_t k = 0; k < mol.nAtoms; k++) {
        const Atom &A_atom = mol.atoms[k];
        double ZA = static_cast<double>(mol.atomicNumbers(k));
        if (ZA > 2){
            ZA -= 2; // counting only the valence
        }
        double gammaAA = gamma(k, k);

        for (const auto &ao_A : A_atom.mAOs) {
            string aoType = string(1, ao_A.AO_type.back());
            double IA = params[A_atom.name].ionizationAffinity[aoType];

            // H_core

            double sumZBgammaAB = 0.0;
            size_t j_AO = 0;
            for (size_t j = 0; j < mol.nAtoms; j++) {
                const Atom &B_atom = mol.atoms[j];
                double ZB = static_cast<double>(mol.atomicNumbers(j));
                if (ZB > 2){
                    ZB -= 2; // counting only the valence
                }

                double beta_A = params[A_atom.name].beta;
                double beta_B = params[B_atom.name].beta;
                double beta_avg = -(beta_A + beta_B) / 2.0;

                if (k != j) {
                    sumZBgammaAB += ZB * gamma(k, j);
                }

                for (const auto &ao_B : B_atom.mAOs) {
                    if (k_AO != j_AO) {
                        H_core(k_AO, j_AO) = beta_avg * S(k_AO, j_AO);
                    }
                    j_AO++;
                }
            }

            H_core(k_AO, k_AO) = -(IA) - (ZA - 0.5) * gammaAA - sumZBgammaAB;
            k_AO++;
        }
    }
}

void CNDO::calculateFockMatrix() {
    size_t u_AO = 0; // Index for atomic orbitals
    for (size_t A = 0; A < mol.nAtoms; ++A) {
        Atom &atom_A = mol.atoms[A];
        string elementType = atom_A.name; // Element type of atom A
        double ZA = static_cast<double>(mol.atomicNumbers(A)); // Atomic number of atom A

        for (auto &ao_A : atom_A.mAOs) {
            string aoType = string(1, ao_A.AO_type.back()); // Get the last character of AO_type
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
            for (size_t B = 0; B < mol.nAtoms; ++B) {
                if (B != A) {
                    double gamma_AB = gamma(A, B); // γ_AB value
                    for (size_t v_AO = 0; v_AO < dim; ++v_AO) {
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
            break;
        }
    }

    // Calculate total energy after convergence
    Etotal = calculateTotalEnergy();
    cout << "Total Energy: " << Etotal << " eV" << endl;
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
    Ga.print("Ga: ");
    Gb.print("Gb: ");
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


