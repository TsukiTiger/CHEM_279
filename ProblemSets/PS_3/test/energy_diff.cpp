#include <iostream>
#include "molecule.h"
#include "eht_matrices.h"
#include <armadillo>

using namespace std;
using namespace arma;

int main () {

    cout << "Finding Bond Energy of H2: " << endl << endl;
    Molecule H2("data/H2.txt");  // Assuming data directory is at the same level as the executable

    mat overlap_mat_H2 = calcOverlapMatrix(H2);
    double energy_H2 = calcEnergy(H2, overlap_mat_H2);
    cout << "Energy of H2: " << energy_H2 << " eV" << endl;

    double energy_H = -13.6;
    double H2_bond_energy = energy_H2 - 2 * energy_H;
    cout << "Bond energy of H2: " << H2_bond_energy << " eV" << endl;

    cout << "----------------------------------------" << endl;
    cout << "Finding Energy Difference: C2H2 + H2 <-> C2H4" << endl << endl;

    Molecule C2H2("data/C2H2.txt");
    Molecule C2H4("data/C2H4.txt");
    Molecule H2_2("data/H2.txt");

    mat overlap_mat_C2H2 = calcOverlapMatrix(C2H2);
    double energy_C2H2 = calcEnergy(C2H2, overlap_mat_C2H2);

    mat overlap_mat_C2H4 = calcOverlapMatrix(C2H4);
    double energy_C2H4 = calcEnergy(C2H4, overlap_mat_C2H4);

    mat overlap_mat_H2_2 = calcOverlapMatrix(H2_2);
    double energy_H2_2 = calcEnergy(H2_2, overlap_mat_H2_2);

    cout << "Energy of C2H4: " << energy_C2H4 << " eV" << endl;
    cout << "Energy of C2H2: " << energy_C2H2 << " eV" << endl;
    cout << "Energy of H2: " << energy_H2_2 << " eV" << endl << endl;

    double C2H4_energy_diff = energy_C2H4 - energy_C2H2 - energy_H2_2;
    C2H4_energy_diff *= 96.485;  // Conversion factor to kJ/mol

    cout << "Energy difference for the chemical reaction C2H2 + H2 <-> C2H4: "
         << C2H4_energy_diff << " kJ/mol" << endl;  // Corrected the unit

    return 0;
}
