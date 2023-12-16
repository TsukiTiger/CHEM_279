//
// Created by Chongye Feng on 12/10/23.
//
#include <iostream>
#include "../inc/Molecule.h"
#include "../inc/CNDO.h"

using namespace std;

const int MAX_ITER = 5000;
const double TOL = 1e-6;

void test_H2() {
    cout << "Testing H2 molecule..." << endl;
    Molecule h2Molecule("data/H2.txt");

    cout << "molecule done..." << endl;

    // Initialize and run CNDO simulation for H2
    CNDO cndoH2(h2Molecule, MAX_ITER, TOL);
    cout << "CNDO/2 done..." << endl;
    cout << "H_core before including repulsion correction:" << endl;
    cndoH2.printResults();

    cndoH2.buildCoreHamiltonian(true);
    cout << "H_core AFTER including repulsion correction:" << endl;

    cndoH2.printResults();

    cndoH2.runSCF();

    // Print results in a format similar to .out file
    cndoH2.printResults();
}

void test_OH() {
    cout << "Testing OH radical..." << endl;
    Molecule ohMolecule("data/OH.txt");

    // Initialize and run CNDO simulation for OH
    CNDO cndoOH(ohMolecule, MAX_ITER, TOL);
    cndoOH.runSCF();

    // Print results in a format similar to .out file
    cndoOH.printResults();
}

void test_C2H4() {
    cout << "Testing C2H4 molecule..." << endl;
    Molecule C2H4Molecule("data/C2H4.txt");

    cout << "molecule done..." << endl;

    // Initialize and run CNDO simulation for H2
    CNDO cndoC2H4(C2H4Molecule, MAX_ITER, TOL);
    cout << "CNDO/2 done..." << endl;
    cout << "H_core before including repulsion correction:" << endl;
    cndoC2H4.runSCF();

    cndoC2H4.buildCoreHamiltonian(true);
    cout << "H_core AFTER including repulsion correction:" << endl;

    cndoC2H4.runSCF();


//    cndoC2H4.printResults();
    double refEnergy = 0.0433634 * -20.1;
    cout << "Ref Energy: " << refEnergy << " eV" << endl;

    // Print results in a format similar to .out file
//    cndoC2H4.printResults();
}

void test_HF() {
    cout << "Testing HF+ molecule..." << endl;
    Molecule HFMolecule("data/HF+.txt");

    cout << "molecule done..." << endl;

    // Initialize and run CNDO simulation for HF+
    CNDO cndoHF(HFMolecule, MAX_ITER, TOL);
    cout << "CNDO/2 done..." << endl;

    cndoHF.runSCF();

    // Print results in a format similar to .out file
    cndoHF.printResults();
}

int main() {
//     test_H2(); // checked Done
    // test_OH();
     test_C2H4();
//    test_HF();
    return 0;
}
