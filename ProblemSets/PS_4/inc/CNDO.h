#include <iostream>
#include <vector>
#include <armadillo>
#include <unordered_map>
#include <string>
#include "molecule.h" // Molecule and BasisFunction class

using namespace std;
using namespace arma;

// Constants
const double HARTREE_TO_EV = 27.211396641308;

// Structure to store CNDO parameters
struct CNDOParams {
    map<string, double> ionizationAffinity; // Ionization potentials and electron affinities
    double beta; // Beta parameter
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
    void buildCoreHamiltonian();

    // Fock matrix calculation
    void calculateFockMatrix();

    // SCF procedure
    void runSCF();

    // Total energy calculation
    double calculateTotalEnergy() ;

    //Initialize the Parameters for the CNDO model
    void initializeCNDOParams();

    void printResults() const;

private:
    Molecule &mol; // Reference to the molecular basis
    int max_iter; // Maximum number of SCF iterations
    double tol; // Convergence tolerance for SCF

    int dim; // Dimension of the matrices (number of atomic orbitals)
    mat Pa, Pb; // Alpha and Beta density matrices
    mat Ga, Gb; // Auxiliary matrices for calculations (if needed)
    mat Ca, Cb; // Coefficient matrices for alpha and beta spins
    vec Ea, Eb; // Eigenvalues (orbital energies) for alpha and beta spins
    mat H_core; // Core Hamiltonian matrix
    mat S; // Overlap matrix
    mat F_alpha, F_beta; // Alpha and Beta Fock matrices
    mat gamma; // Two-electron integral matrix (gamma values)

    int q, p; // Number of alpha and beta electrons
    double Ee, Etotal, Ec; // Electronic, total, and core energies

    unordered_map<string, CNDOParams> params; // CNDO parameters
    void initializeParams();
    void updateG();
};
