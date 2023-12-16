//
// Created by Chongye Feng on 12/6/23.
//

#ifndef EHT_FINAL_UTILS_H
#define EHT_FINAL_UTILS_H


#include "BasisFunction.h"
#include "Molecule.h"
#include <armadillo>

class BasisFunction;

double I2e_pG(const arma::rowvec &Ra, const arma::rowvec &Rb, double sigmaA, double sigmaB);
double eval_2eI_s(const BasisFunction &ao1, const BasisFunction &ao2);
double factorial(int n);
double binomialCoef(int m, int n);
double doubleFactorial(int n);


#endif //EHT_FINAL_UTILS_H
