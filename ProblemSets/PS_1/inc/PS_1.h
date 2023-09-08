// File: /inc/PS_1.h

#pragma once

#include <vector>
#include "../utils/atom.h"
#include "../utils/vector3d.h"

double calculate_distance(const Atom& atom_i, const Atom& atom_j);
double calculate_lj_energy_au(const Atom& atom_i, const Atom& atom_j);
bool read_file(const std::string& path, std::vector<Atom>& atoms);
double calculate_total_energy_au(const std::vector<Atom>& atoms);
std::vector<Vector3D> calculate_analytical_forces(const std::vector<Atom>& atoms);
std::vector<Vector3D> calculate_forward_difference_force(const std::vector<Atom>& atoms, double h);
std::vector<Vector3D> calculate_central_difference_force(const std::vector<Atom>& atoms, double h);
double line_search(const std::vector<Atom>& atoms, const std::vector<Vector3D>& forces);
bool is_converged(double energy_old, double energy_new, const std::vector<Vector3D>& forces);
void optimize_structure(std::vector<Atom>& atoms);
