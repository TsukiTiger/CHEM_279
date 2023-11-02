#include "eht_matrices.h"

// Assuming constant_K is global and static to this file.
// You can update this value or pass it as an argument to the relevant function.
static const double constant_K = 1.75; // Assign the appropriate value for constant_K

mat calcHamiltonianMatrix(const Molecule& molecule, const mat& overlap_mat) {
    // Initialize the Hamiltonian matrix
    mat hamiltonian_mat = arma::zeros(molecule.nBasisFunctions, molecule.nBasisFunctions);

    // Define a map to store the atom type to Hamiltonian value mappings
    map<string, double> hamiltonianValues = {
            {"H-1s", -13.6},
            {"C-2s", -21.4},
            {"C-2px", -11.4},
            {"C-2py", -11.4},
            {"C-2pz", -11.4}
    };

    // Loop over all the basis functions of the molecule
    for (int i = 0; i < molecule.nBasisFunctions; i++) {
        for (int j = 0; j < molecule.nBasisFunctions; j++) {
            string atomType1 = molecule.basisFunctionsList[i].AO_type;
            string atomType2 = molecule.basisFunctionsList[j].AO_type;

            // Get the Hamiltonian values given the atom types (i.e H-1s, C-2s, etc.)
            double h_i = hamiltonianValues[atomType1];
            double h_j = hamiltonianValues[atomType2];

            // Assign diagonal elements
            if (i == j) {
                hamiltonian_mat(i, j) = h_i;
            } else {
                // Calculate off-diagonal elements
                hamiltonian_mat(i, j) = constant_K / 2 * (h_i + h_j) * overlap_mat(i, j);
            }
        }
    }
    return hamiltonian_mat;
}

mat calcXMatrix(const mat& overlap_mat) {
    // Initialize the X matrix
    mat X_mat = arma::zeros(overlap_mat.n_rows, overlap_mat.n_cols);

    // Diagonalize the overlap matrix
    vec s_eigval;
    mat U_eigvec;
    eig_sym(s_eigval, U_eigvec, overlap_mat);

    // Get inverse of square root of eigenvalues
    vec s_eigval_inv_sqrt = 1.0 / sqrt(s_eigval);

    // Diagonalize s_eigval_inv_sqrt
    mat s_diag_mat = diagmat(s_eigval_inv_sqrt);

    // Calculate X matrix
    X_mat = U_eigvec * s_diag_mat * U_eigvec.t();

    return X_mat;
}

mat calcHamiltonianPrimeMatrix(const mat& X_mat, const mat& hamiltonian_mat) {
    // Initialize the Hamiltonian Prime matrix
    mat h_prime_mat = arma::zeros(hamiltonian_mat.n_rows, hamiltonian_mat.n_cols);

    // Calculate Hamiltonian Prime matrix
    h_prime_mat = X_mat.t() * hamiltonian_mat * X_mat;

    return h_prime_mat;
}

double calcEnergy(const mat& X_mat, const mat& hamiltonian_prime_mat, int nElectrons) {
    // Diagonalize the Hamiltonian Prime matrix
    vec energy_eigval;
    mat V_eigvec;
    eig_sym(energy_eigval, V_eigvec, hamiltonian_prime_mat);

    // Form the MO coffecient matrix
    mat C_mat = X_mat * V_eigvec;

    // Calculate the energy by summing over the occupied orbitals
    double energy = 0;
    for (int i = 0; i < nElectrons; i++) {
        energy += 2 * energy_eigval(i);
    }

    return energy;
}

double calcEnergy(const Molecule& molecule, const mat& overlap_mat) {
    mat hamiltonian_mat = calcHamiltonianMatrix(molecule, overlap_mat);
    mat X_mat = calcXMatrix(overlap_mat);
    mat hamiltonian_prime_mat = calcHamiltonianPrimeMatrix(X_mat, hamiltonian_mat);
    double energy = calcEnergy(X_mat, hamiltonian_prime_mat, molecule.nElectrons);

    return energy;
}
