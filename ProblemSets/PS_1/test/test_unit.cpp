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

void test_cases() {
    std::cout << "Testing custom cases..." << std::endl;
    std::cout << "---------------------------------------" << std::endl;

    // Test Case 1
    std::vector<Atom> atoms1;
    read_file("data/test_case_1.txt", atoms1);
    double energy1 = calculate_total_energy_au(atoms1);
    std::cout << "Test Case 1 - Calculated Energy: " << energy1 << std::endl;
    std::cout << "Test Case 1 - Reference Energy: -1.7358" << std::endl;

    auto forces1 = calculate_analytical_forces(atoms1);
    std::vector<Vector3D> ref_forces1 = {
            {-0.042191, 0.56398, 0},
            {-0.89489, 0, 0},
            {-0.042191, -0.56398, 0},
            {0.97927, 0, 0}
    };

    double sum_of_square_forces = 0.0;

    for (size_t i = 0; i < forces1.size(); ++i) {
        sum_of_square_forces += forces1[i].x * forces1[i].x + forces1[i].y * forces1[i].y + forces1[i].z * forces1[i].z;
        std::cout << "Atom " << i + 1 << " Force: (" << forces1[i].x << ", " << forces1[i].y << ", " << forces1[i].z << ")" << std::endl;
        std::cout << "Atom " << i + 1 << " Reference Force: (" << ref_forces1[i].x << ", " << ref_forces1[i].y << ", " << ref_forces1[i].z << ")" << std::endl;
    }
    double force_norm1 = std::sqrt(sum_of_square_forces);
    std::cout << "Force norm for Test Case 1: " << force_norm1 << " Reference: 1.549" << std::endl;
    std::cout << "---------------------------------------" << std::endl;

    // Test Case 2
    std::vector<Atom> atoms2;
    read_file("data/test_case_2.txt", atoms2);
    double energy2 = calculate_total_energy_au(atoms2);
    std::cout << "Test Case 2 - Calculated Energy: " << energy2 << std::endl;
    std::cout << "Test Case 2 - Reference Energy: -2.0192" << std::endl;

    auto forces2 = calculate_analytical_forces(atoms2);
    std::vector<Vector3D> ref_forces2 = {
            {0, 2.155, 0},
            {0, -1.6323, 0},
            {0, -0.52269, 0}
    };

    sum_of_square_forces = 0.0;
    for (size_t i = 0; i < forces2.size(); ++i) {
        sum_of_square_forces += forces2[i].x * forces2[i].x + forces2[i].y * forces2[i].y + forces2[i].z * forces2[i].z;
        std::cout << "Atom " << i + 1 << " Force: (" << forces2[i].x << ", " << forces2[i].y << ", " << forces2[i].z << ")" << std::endl;
        std::cout << "Atom " << i + 1 << " Reference Force: (" << ref_forces2[i].x << ", " << ref_forces2[i].y << ", " << ref_forces2[i].z << ")" << std::endl;
    }
    double force_norm2 = std::sqrt(sum_of_square_forces);;
    std::cout << "Force norm for Test Case 2: " << force_norm2 << " Reference: 2.7534" << std::endl;
    std::cout << "---------------------------------------" << std::endl;
}

int main() {
    test_read_file();
    test_simple_systems();
    test_structure_optimization();
    test_force_approximations();
    plot_error_vs_step_size();
    test_cases();
    return 0;
}

