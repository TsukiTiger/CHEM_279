#include <iostream>
#include "molecule.h"
#include "eht_matrices.h"

using namespace std;  // Adds std namespace for cout, endl, etc.

int main(int argc, char* argv[]) {
    if (argc != 2) {
        cerr << "Usage: " << argv[0] << " <filename>\n";
        return 1;
    }

    string filename = argv[1];  // Convert to std::string if Molecule constructor requires it

    Molecule molecule(filename);

    molecule.printMoleculeInfo();

    cout << "----------------------------------------" << endl << endl;

    // Test contracted overlap integral
    mat overlap_mat = calcOverlapMatrix(molecule);
    cout << "Overlap matrix:\n" << overlap_mat << endl;

    // Test Hamiltonian matrix
    mat hamiltonian_mat = calcHamiltonianMatrix(molecule, overlap_mat);
    cout << "Hamiltonian matrix:\n" << hamiltonian_mat << endl;

    // Test X matrix
    mat X_mat = calcXMatrix(overlap_mat);
    cout << "X matrix:\n" << X_mat << endl;

    // Test Hamiltonian Prime matrix
    mat hamiltonian_prime_mat = calcHamiltonianPrimeMatrix(X_mat, hamiltonian_mat);
    cout << "Hamiltonian Prime matrix:\n" << hamiltonian_prime_mat << endl;

    // Test calcEnergy
    double energy = calcEnergy(molecule, overlap_mat);
    cout << "Energy: " << energy << endl;

    return 0;

    return 0;
}
