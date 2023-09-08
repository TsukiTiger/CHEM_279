#include <iostream>
#include <vector>
#include "../inc/PS_1.h"

void test_read_file() {
    std::vector<Atom> atoms;

    // Test the valid file first
    if (read_file("data/coords.txt", atoms)) {
        std::cout << "Test with coords.txt passed!" << std::endl;
    } else {
        std::cout << "Unable to read data from coords.txt!" << std::endl;
    }

    // Test the erroneous file
    if (read_file("data/coords_err.txt", atoms)) {
        std::cout << "Test with coords_err.txt passed unexpectedly!" << std::endl;
    } else {
        std::cout << "Test with coords_err.txt failed as expected!" << std::endl;
    }

    std::cout << "---------------------------------------" << std::endl;
}

void test_simple_systems() {
    Atom atom1 = {79,0.0, 0.0, 0.0}; // Assuming your Atom structure can be initialized like this
    Atom atom2 = {79,2.0, 0.0, 0.0}; // Separation of 2 units

    double energy = calculate_lj_energy_au(atom1, atom2);
    std::cout << "Energy for Au2 separated by 2 units: " << energy << std::endl;
    // You should compare this with the known value or reference.

    std::cout << "---------------------------------------" << std::endl;
}

void test_structure_optimization() {
    std::vector<Atom> atoms = {
            {79, 0.0, 0.0, 0.0},
            {79, 2.5, 0.0, 0.0} // Just an initial configuration
    };

    std::cout << "Energy before optimizing: " << calculate_total_energy_au(atoms) << std::endl;

    optimize_structure(atoms);

    // Print out the optimized coordinates and energy:
    for (const auto& atom : atoms) {
        std::cout << "Atom coordinates: (" << atom.x << ", " << atom.y << ", " << atom.z << ")" << std::endl;
    }
    std::cout << "Optimized energy: " << calculate_total_energy_au(atoms) << std::endl;

    std::cout << "---------------------------------------" << std::endl;
}

void test_force_approximations() {
    std::vector<Atom> atoms = {
            {79, 0.0, 0.0, 0.0},
            {79, 2.5, 0.0, 0.0}
    };

    // Here, assuming you have functions like calculate_forward_difference_force and
    // calculate_central_difference_force and calculate_analytical_forces.
    Vector3D analytical_force = calculate_analytical_forces(atoms)[0];
    Vector3D forward_diff_force = calculate_forward_difference_force(atoms, 0.01)[0];
    Vector3D central_diff_force = calculate_central_difference_force(atoms, 0.01)[0];

    std::cout << "Analytical force on atom 1: (" << analytical_force.x << ", "
              << analytical_force.y << ", " << analytical_force.z << ")" << std::endl;
    std::cout << "Forward difference approximated force on atom 1: (" << forward_diff_force.x << ", "
              << forward_diff_force.y << ", " << forward_diff_force.z << ")" << std::endl;
    std::cout << "Central difference approximated force on atom 1: (" << central_diff_force.x << ", "
              << central_diff_force.y << ", " << central_diff_force.z << ")" << std::endl;

    // Similarly, compare the magnitudes, and ensure that they are close to each other.
    std::cout << "---------------------------------------" << std::endl;
}

void plot_error_vs_step_size() {
    std::vector<Atom> atoms = {
            {79, 0.0, 0.0, 0.0},
            {79, 2.5, 0.0, 0.0}
    };

    std::vector<double> step_sizes = {0.1, 0.01, 0.001, 0.0001};

    for (double h : step_sizes) {
        Vector3D analytical_force = calculate_analytical_forces(atoms)[0];
        Vector3D forward_diff_force = calculate_forward_difference_force(atoms, h)[0];
        Vector3D central_diff_force = calculate_central_difference_force(atoms, h)[0];

        double forward_error = std::sqrt(std::pow(analytical_force.x - forward_diff_force.x, 2) +
                                         std::pow(analytical_force.y - forward_diff_force.y, 2) +
                                         std::pow(analytical_force.z - forward_diff_force.z, 2));

        double central_error = std::sqrt(std::pow(analytical_force.x - central_diff_force.x, 2) +
                                         std::pow(analytical_force.y - central_diff_force.y, 2) +
                                         std::pow(analytical_force.z - central_diff_force.z, 2));


        std::cout << "h: " << h << ", Forward difference error: " << forward_error
                  << ", Central difference error: " << central_error << std::endl;
    }

    // Typically, you would use this data to create a log-log plot and examine the slope.
    std::cout << "---------------------------------------" << std::endl;
}


int main() {
    test_read_file();
    test_simple_systems();
    test_structure_optimization();
    test_force_approximations();
    plot_error_vs_step_size();
    return 0;
}

