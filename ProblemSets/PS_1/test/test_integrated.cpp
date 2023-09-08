//
// Created by Chongye Feng on 9/6/23.
//
#include <iostream>
#include "../inc/PS_1.h"

int main() {
    std::vector<Atom> atoms;

    if (!read_file("data/coords.txt", atoms)) {
        std::cerr << "Error reading coordinates file." << std::endl;
        return 1;
    }

    std::cout << "Atomic positions before optimization:" << std::endl;
    for (const Atom& atom : atoms) {
        std::cout << "Atom " << atom.atomic_number << " at ("
                  << atom.x << ", " << atom.y << ", " << atom.z << ")" << std::endl;
    }

    std::cout << "Energy before optimizing: " << calculate_total_energy_au(atoms) << std::endl;
    std::cout << "---------------------------------------" << std::endl;

    optimize_structure(atoms);

    std::cout << "Optimized atomic positions:" << std::endl;
    for (const Atom& atom : atoms) {
        std::cout << "Atom " << atom.atomic_number << " at ("
                  << atom.x << ", " << atom.y << ", " << atom.z << ")" << std::endl;
    }

    std::cout << "Energy after optimizing: " << calculate_total_energy_au(atoms) << std::endl;
    std::cout << "---------------------------------------" << std::endl;

    return 0;
}
