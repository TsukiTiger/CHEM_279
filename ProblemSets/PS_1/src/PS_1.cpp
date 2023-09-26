#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include "../utils/constants.h"
#include "../utils/atom.h"
#include "../utils/vector3d.h"
#include "../inc/PS_1.h"  // Include the new header file

double epsilon = std::sqrt(constant::epsilon_au * constant::epsilon_au);
double sigma = std::sqrt(constant::sigma_au * constant::sigma_au);
const double energy_tol = 1e-2;
const double force_tol = 1e-6;


double calculate_distance(const Atom& atom_i, const Atom& atom_j) {
    return std::sqrt((atom_i.x - atom_j.x) * (atom_i.x - atom_j.x) +
                     (atom_i.y - atom_j.y) * (atom_i.y - atom_j.y) +
                     (atom_i.z - atom_j.z) * (atom_i.z - atom_j.z));
}

double calculate_lj_energy_au(const Atom& atom_i, const Atom& atom_j) {
    double rij = calculate_distance(atom_i, atom_j);
    double term1 = std::pow(sigma / rij, 12);
    double term2 = std::pow(sigma / rij, 6);
    return epsilon * (term1 - 2 * term2);
}

// Function to read atomic coordinates from a file
bool read_file(const std::string& path, std::vector<Atom>& atoms) {
    atoms.clear();
    Atom atom;
    int num_atoms;

    std::ifstream input_file(path);
    if (!input_file) {
        std::cerr << "Error: Could not open input file." << std::endl;
        return false;
    }

    input_file >> num_atoms;

    while (input_file >> atom.atomic_number >> atom.x >> atom.y >> atom.z) {
        if (atom.atomic_number != 79) {
            std::cerr << "Error: Only Au atoms are supported." << std::endl;
            return false;
        }
        atoms.push_back(atom);
    }

    if (atoms.size() != num_atoms) {
        std::cerr << "Error: Mismatch in declared and actual number of atoms." << std::endl;
        return false;
    }

    input_file.close();
    return true;
}

// Function to calculate the total LJ potential energy for a system of atoms
double calculate_total_energy_au(const std::vector<Atom>& atoms) {
    double total_energy = 0.0;
    for (size_t i = 0; i < atoms.size(); ++i) {
        for (size_t j = i + 1; j < atoms.size(); ++j) {
            total_energy += calculate_lj_energy_au(atoms[i], atoms[j]);
        }
    }

    return total_energy;
}

// Function to calculate analytical LJ forces
std::vector<Vector3D> calculate_analytical_forces(const std::vector<Atom>& atoms) {
    std::vector<Vector3D> forces(atoms.size(), Vector3D(0.0, 0.0, 0.0));

    for (size_t i = 0; i < atoms.size(); ++i) {
        for (size_t j = 0; j < atoms.size(); ++j) {
            if (i != j) {
                const Atom& atom_i = atoms[i];
                const Atom& atom_j = atoms[j];

                double rij = calculate_distance(atom_i, atom_j);
                double r13 = std::pow(rij, 13);
                double r7 = std::pow(rij, 7);

                // Calculating the magnitude of the force
                double force_mag = epsilon * (12 * std::pow(sigma, 12) / r13 - 12 * std::pow(sigma, 6) / r7);

                // Calculating the direction of the force
                double dx = atom_i.x - atom_j.x;
                double dy = atom_i.y - atom_j.y;
                double dz = atom_i.z - atom_j.z;

                // Updating the forces for atom i
                forces[i].x += force_mag * dx / rij;
                forces[i].y += force_mag * dy / rij;
                forces[i].z += force_mag * dz / rij;
            }
        }
    }
    return forces;
}

std::vector<Vector3D> calculate_forward_difference_force(const std::vector<Atom>& atoms, double h) {
    std::vector<Vector3D> forces(atoms.size());

    for (size_t i = 0; i < atoms.size(); ++i) {
        std::vector<Atom> perturbed_atoms = atoms;

        perturbed_atoms[i].x += h;
        double energy_plus_hx = calculate_total_energy_au(perturbed_atoms);

        perturbed_atoms[i].x -= h; // Reset the perturbation
        perturbed_atoms[i].y += h;
        double energy_plus_hy = calculate_total_energy_au(perturbed_atoms);

        perturbed_atoms[i].y -= h; // Reset the perturbation
        perturbed_atoms[i].z += h;
        double energy_plus_hz = calculate_total_energy_au(perturbed_atoms);
        double energy_initial = calculate_total_energy_au(atoms);

        forces[i].x = -(energy_plus_hx - energy_initial) / h;
        forces[i].y = -(energy_plus_hy - energy_initial) / h;
        forces[i].z = -(energy_plus_hz - energy_initial) / h;
    }

    return forces;
}

std::vector<Vector3D> calculate_central_difference_force(const std::vector<Atom>& atoms, double h) {
    std::vector<Vector3D> forces(atoms.size());

    for (size_t i = 0; i < atoms.size(); ++i) {
        std::vector<Atom> perturbed_atoms_plus = atoms;
        std::vector<Atom> perturbed_atoms_minus = atoms;

        perturbed_atoms_plus[i].x += h;
        perturbed_atoms_minus[i].x -= h;
        double energy_plus_hx = calculate_total_energy_au(perturbed_atoms_plus);
        double energy_minus_hx = calculate_total_energy_au(perturbed_atoms_minus);

        perturbed_atoms_plus[i].x -= h; // Reset the perturbation
        perturbed_atoms_plus[i].y += h;
        perturbed_atoms_minus[i].y -= h;
        double energy_plus_hy = calculate_total_energy_au(perturbed_atoms_plus);
        double energy_minus_hy = calculate_total_energy_au(perturbed_atoms_minus);

        perturbed_atoms_plus[i].y -= h; // Reset the perturbation
        perturbed_atoms_plus[i].z += h;
        perturbed_atoms_minus[i].z -= h;
        double energy_plus_hz = calculate_total_energy_au(perturbed_atoms_plus);
        double energy_minus_hz = calculate_total_energy_au(perturbed_atoms_minus);

        forces[i].x = -(energy_plus_hx - energy_minus_hx) / (2 * h);
        forces[i].y = -(energy_plus_hy - energy_minus_hy) / (2 * h);
        forces[i].z = -(energy_plus_hz - energy_minus_hz) / (2 * h);
    }

    return forces;
}

// Function to perform a 1D line search to find the optimal step size
double line_search(const std::vector<Atom>& atoms, const std::vector<Vector3D>& forces) {
    double alpha = 0.1; // Initial step size (PS: I tried 0.01, it would be too small for it to converge)
    double energy_initial = calculate_total_energy_au(atoms);
    double energy_new;
    const double c = 0.5; // Line search parameter

    while (alpha > 1e-6) { // Convergence based on step size
        // Make a trial move in the steepest descent direction
        std::vector<Atom> trial_atoms = atoms;
        for (size_t i = 0; i < atoms.size(); ++i) {
            trial_atoms[i].x -= alpha * forces[i].x;
            trial_atoms[i].y -= alpha * forces[i].y;
            trial_atoms[i].z -= alpha * forces[i].z;
        }

        // Calculate the energy at the trial position
        energy_new = calculate_total_energy_au(trial_atoms);

        // Check if the energy decreased sufficiently
        if (energy_new < energy_initial - alpha * c) {
            return alpha; // Found a suitable step size
        }

        alpha *= 0.5; // Reduce the step size
    }

    return alpha; // Return the last step size
}

// Check convergence based on energy change and max force component
bool is_converged(double energy_old, double energy_new, const std::vector<Vector3D>& forces) {
    double max_force = 0.0;
    for (const auto& force : forces) {
        max_force = std::max(max_force, force.norm()); // Compute the magnitude of each force
    }
    return (std::abs(energy_new - energy_old) < energy_tol) || (max_force < force_tol);
}

void optimize_structure(std::vector<Atom>& atoms) {
    const int max_iterations = 5000;
    double energy_old = calculate_total_energy_au(atoms);

    for (int iter = 0; iter < max_iterations; ++iter) {
        // Calculate forces and energy
        std::vector<Vector3D> forces = calculate_analytical_forces(atoms);

        // Perform a line search to find the optimal step size
        double alpha = line_search(atoms, forces);

        // Update atomic positions
        for (size_t i = 0; i < atoms.size(); ++i) {
            atoms[i].x -= alpha * forces[i].x;
            atoms[i].y -= alpha * forces[i].y;
            atoms[i].z -= alpha * forces[i].z;
        }

        // Calculate the new energy
        double energy_new = calculate_total_energy_au(atoms);

        // Debug outputs
        std::cout << "Iteration " << iter << ": Energy = " << energy_new
                  << ", Max Force = " << forces[0].norm()
                  << ", Step Size = " << alpha << std::endl;

        // Check for convergence
        if (is_converged(energy_old, energy_new, forces)) {
            std::cout << "Converged after " << iter << " iterations." << std::endl;
            return;
        }

        energy_old = energy_new;
    }

    std::cerr << "Optimization did not converge within the maximum allowed iterations!" << std::endl;
}
