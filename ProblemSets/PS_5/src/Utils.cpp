//
// Created by Chongye Feng on 12/6/23.
//

#include "../inc/Utils.h"

double I2e_pG(const arma::rowvec &Ra, const arma::rowvec &Rb, double sigmaA, double sigmaB)
{
    double U = pow(M_PI*sigmaA, 1.5) * pow(M_PI*sigmaB, 1.5); //eq 3.8 and 3.11
    double V2= 1.0/(sigmaA+sigmaB); //eq 3.9
    double Rd = arma::norm(Ra-Rb, 2);
    if (Rd==0.0) {
        //use eq 3.15
        return U * 2 * sqrt(V2/M_PI);
    }
    //eq 3.14 need sqrt T
    double srT = sqrt(V2) * Rd;
    //eq 3.14
    double result = U / Rd * erf(srT);
    return result;
}

double eval_2eI_s(const BasisFunction &ao1, const BasisFunction &ao2) {
    // Check if both ao1 and ao2 are s orbitals (lmn vectors should be {0,0,0})

    if (!(arma::accu(ao1.getLmn()) == 0) or !(arma::accu(ao2.getLmn()) == 0)) {
        std::cerr << "Error: Only s orbitals are allowed for eval_2eI_s function." << std::endl;
        return 0.0;
    }

    // Get the contraction coefficients multiplied by the normalization constants
    // Assuming da and db are arma::vec
    arma::vec da(ao1.getContractionCoeffs().size());
    arma::vec db(ao2.getContractionCoeffs().size());

    for (size_t i = 0; i < da.size(); ++i) {
        da[i] = ao1.getContractionCoeffs()[i] * ao1.getNormConstants()[i];
        db[i] = ao2.getContractionCoeffs()[i] * ao2.getNormConstants()[i];
    }


    // Extract exponents and centers
    const arma::vec& alphaa = ao1.getExponents();
    const arma::vec& alphab = ao2.getExponents();
    const arma::rowvec& Ra = ao1.getCenter();
    const arma::rowvec& Rb = ao2.getCenter();

    // Assuming the length of alpha vectors is 3 for STO-3G basis set
    const int len = 3;

    double sum = 0.0;
    // Nested loops to compute the sum over all primitives in the contracted Gaussians
    for (int k1 = 0; k1 < len; ++k1) {
        for (int k2 = 0; k2 < len; ++k2) {
            double sigmaA = 1.0 / (alphaa(k1) + alphaa(k2)); // Pre-factor for AO1

            for (int j1 = 0; j1 < len; ++j1) {
                for (int j2 = 0; j2 < len; ++j2) {
                    double sigmaB = 1.0 / (alphab(j1) + alphab(j2)); // Pre-factor for AO2
                    double I2e = I2e_pG(Ra, Rb, sigmaA, sigmaB); // Two-electron integral calculation

                    // Accumulate the result with contraction and normalization factors
                    sum += da(k1) * da(k2) * db(j1) * db(j2) * I2e;
                }
            }
        }
    }

    return sum; // Return the computed two-electron integral (gamma value)
}

double factorial(int n) {
    int result = 1;
    for (int i = n; i > 0; i--) {
        result *= i;
    }
    return result;
}

double binomialCoef(int m, int n) {
    if (m < n) {
        return 0;
    }
    return factorial(m) / (factorial(n) * factorial(m - n));
}

double doubleFactorial(int n) {
    int result = 1;
    for (int i = n; i > 0; i -= 2) {
        result *= i;
    }
    return result;
}