//
// Created by Chongye Feng on 12/6/23.
//

#ifndef BASISFUNCTION_H
#define BASISFUNCTION_H

#include <armadillo>
#include <vector>
#include <string>
#include "Utils.h"
#include "Molecule.h"

class BasisFunction {
public:
    BasisFunction(const std::string& AO_type,
                  const arma::rowvec& center,
                  const arma::vec& lmn,
                  const arma::vec& exponents,
                  const arma::vec& contraction_coeffs);

    const std::string &getAoType() const;
    const arma::rowvec &getCenter() const;
    const arma::vec &getLmn() const;
    const arma::vec &getExponents() const;
    const arma::vec &getContractionCoeffs() const;
    const std::vector<double> &getNormConstants() const;

    void setAoType(const std::string &aoType);
    void setCenter(const arma::rowvec &center);
    void setLmn(const arma::vec &lmn);
    void setExponents(const arma::vec &exponents);
    void setContractionCoeffs(const arma::vec &contractionCoeffs);
    void setNormConstants(const std::vector<double> &normConstants);

private:
    void calcNormConstants();

    // Member variables
    std::string AO_type;
    arma::rowvec center;
    arma::vec lmn;
    arma::vec exponents;
    arma::vec contraction_coeffs;
    std::vector<double> norm_constants;

};

#endif // BASISFUNCTION_H
