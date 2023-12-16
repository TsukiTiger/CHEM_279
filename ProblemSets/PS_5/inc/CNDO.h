//
// Created by Chongye Feng on 12/6/23.
//

#ifndef EHT_FINAL_CNDO_H
#define EHT_FINAL_CNDO_H


#include <iostream>
#include <vector>
#include <armadillo>
#include <unordered_map>
#include <string>
#include "Molecule.h"
#include "Utils.h"

using namespace std;
using namespace arma;

// Constants
const double HARTREE_TO_EV = 27.211396641308;

// Structure to store CNDO parameters
struct CNDOParams {
    map<string, double> ionizationAffinity; // Ionization potentials and electron affinities
    double beta; // Beta parameter
};

struct RepulsionParams {
    double alpha;
    double gamma;
    double omega;
    double r;
};

// Function to initialize CNDO parameters
unordered_map<string, CNDOParams> initializeCNDOParams();

// Define the CNDO class for calculations
class CNDO {
public:
    // Constructor
    CNDO(Molecule &mol_i, int max_it, double tolerance);

    // Initialization and setup
    void initialize();

    // Core Hamiltonian calculation
    void buildCoreHamiltonian(bool repulsionCorrection);

    // Fock matrix calculation
    void calculateFockMatrix();

    // SCF procedure
    void runSCF();

    // Total energy calculation
    double calculateTotalEnergy() ;

    //Initialize the Parameters for the CNDO model
    void initializeCNDOParams();
    void rParamsInit();

    void printResults() const;

    void calculateGradient();
    mat calculateNuclearRepulsionGradient();
    mat calculateElectronicGradient();
    void assembleXMatrix();
    void assembleYMatrix();
    void calculateOV_RA();
    void calculateGamma_RA();
    void printGradients() const;

private:
    Molecule &mol; // Reference to the molecular basis
    int max_iter; // Maximum number of SCF iterations
    double tol; // Convergence tolerance for SCF
    bool scfConverged;

    int dim; // Dimension of the matrices (number of atomic orbitals)
    mat Pa, Pb; // Alpha and Beta density matrices
    mat Ga, Gb; // Auxiliary matrices for calculations (if needed)
    mat Ca, Cb; // Coefficient matrices for alpha and beta spins
    vec Ea, Eb; // Eigenvalues (orbital energies) for alpha and beta spins
    mat H_core; // Core Hamiltonian matrix
    mat S; // Overlap matrix
    mat F_alpha, F_beta; // Alpha and Beta Fock matrices
    mat gamma; // Two-electron integral matrix (gamma values)

    mat gradient;  // Gradient matrix
    mat n_gradient;  // Nuclear Gradient matrix
    mat e_gradient;  // Electronic Gradient matrix
    arma::field<arma::vec> dS_RA;
    arma::field<arma::vec> dGamma_RA;
    mat xMatrix, yMatrix;
    vec OV_vec(const BasisFunction &a, const BasisFunction &b);
    vec d00_vec(const BasisFunction &a, const BasisFunction &b, double sigma_a, double sigma_b);
    vec Gamma_vec(const BasisFunction &a, const BasisFunction &b);
    double d_s_AB(double exp_a, double exp_b, double center_a, double center_b, int lA, int lB);


    int q, p; // Number of alpha and beta electrons
    double Ee, Etotal, Ec; // Electronic, total, and core energies

    unordered_map<string, CNDOParams> params; // CNDO parameters
    unordered_map<std::string, RepulsionParams> rparams;
    double calculateRepulsionCorrection(double RAB, const std::string& atomPair);
    void updateG();
    void updateRepulsionHamiltonian();
};



#endif //EHT_FINAL_CNDO_H
