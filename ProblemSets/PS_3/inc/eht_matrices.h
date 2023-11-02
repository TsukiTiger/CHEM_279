#ifndef EHT_MATRICES_H
#define EHT_MATRICES_H

#include <iostream>
#include <vector>
#include <armadillo>
#include <cmath>
#include <map>
#include "molecule.h"
using namespace arma;
using namespace std;

// Calculate the Hamiltonian matrix from the overlap matrix
mat calcHamiltonianMatrix(const Molecule& molecule, const mat& overlap_mat);

// Calculate the X transformation matrix by diagonalizing the overlap matrix
mat calcXMatrix(const mat& overlap_mat);

// Calculate Hamiltonian Prime matrix from the Hamiltonian matrix and the X matrix
mat calcHamiltonianPrimeMatrix(const mat& X_mat, const mat& hamiltonian_mat);

// Calculate the energy of the molecule given the Hamiltonian Prime matrix
double calcEnergy(const mat& X_mat, const mat& hamiltonian_prime_mat, int nElectrons);

// Calculate the energy of the molecule given the molecule object and the overlap matrix
double calcEnergy(const Molecule& molecule, const mat& overlap_mat);

#endif // EHT_MATRICES_H
