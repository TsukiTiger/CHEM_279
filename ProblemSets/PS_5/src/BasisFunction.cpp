//
// Created by Chongye Feng on 12/6/23.
//

#include "../inc/BasisFunction.h"

#include <armadillo>
#include <vector>
#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <cassert>

BasisFunction::BasisFunction(const std::string& AO_type,
                             const arma::rowvec& center,
                             const arma::vec& lmn,
                             const arma::vec& exponents,
                             const arma::vec& contraction_coeffs)
        : AO_type(AO_type), center(center), lmn(lmn),
          exponents(exponents), contraction_coeffs(contraction_coeffs) {
    calcNormConstants();
}

void BasisFunction::calcNormConstants() {
    // Calculate the overlap integrals for each exponent with itself

    double selfOverlap_1 = Molecule::overlapIntegral3D(center, center, exponents(0), exponents(0), lmn, lmn);
    double selfOverlap_2 = Molecule::overlapIntegral3D(center, center, exponents(1), exponents(1), lmn, lmn);
    double selfOverlap_3 = Molecule::overlapIntegral3D(center, center, exponents(2), exponents(2), lmn, lmn);

    // Calculate the normalization constants
    norm_constants = {1.0 / sqrt(selfOverlap_1), 1.0 / sqrt(selfOverlap_2), 1.0 / sqrt(selfOverlap_3)};
}

const std::string &BasisFunction::getAoType() const {
    return AO_type;
}

const arma::rowvec &BasisFunction::getCenter() const {
    return center;
}

const arma::vec &BasisFunction::getLmn() const {
    return lmn;
}

const arma::vec &BasisFunction::getExponents() const {
    return exponents;
}

const arma::vec &BasisFunction::getContractionCoeffs() const {
    return contraction_coeffs;
}

const std::vector<double> &BasisFunction::getNormConstants() const {
    return norm_constants;
}

void BasisFunction::setAoType(const std::string &aoType) {
    AO_type = aoType;
}

void BasisFunction::setCenter(const arma::rowvec &center) {
    BasisFunction::center = center;
}

void BasisFunction::setLmn(const arma::vec &lmn) {
    BasisFunction::lmn = lmn;
}

void BasisFunction::setExponents(const arma::vec &exponents) {
    BasisFunction::exponents = exponents;
}

void BasisFunction::setContractionCoeffs(const arma::vec &contractionCoeffs) {
    contraction_coeffs = contractionCoeffs;
}

void BasisFunction::setNormConstants(const std::vector<double> &normConstants) {
    norm_constants = normConstants;
}
