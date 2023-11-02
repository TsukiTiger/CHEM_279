#ifndef PS_2_H
#define PS_2_H

#include <cmath>

/**
 * @class GaussianShell
 * @brief Represents a Gaussian shell and contains methods to perform various computations on it.
 */
class GaussianShell {
public:
    double center; ///< Center of the Gaussian shell.
    double exponent; ///< Exponent of the Gaussian shell.
    int angularMomentum; ///< Angular momentum of the Gaussian shell.

    /**
     * @brief Default constructor that initializes a GaussianShell object.
     */
    GaussianShell();


    /*
     * 3 Gaussian object
     * one for each coord (x, y, z)
     *
     * 3D == looping through each coord (x, y, z)
     *
     *
     * */
    /**
     * @brief Parametrized constructor that initializes a GaussianShell object with specific values.
     * @param c Center of the Gaussian shell.
     * @param e Exponent of the Gaussian shell.
     * @param l Angular momentum of the Gaussian shell.
     */
    GaussianShell(double c, double e, int l);

    double computeAnalyticalOverlap_1d(const GaussianShell& shellB) const; ///< Method to compute 1d analytical overlap.
    double computeNumericalOverlap_1d(const GaussianShell& shellB, double lower_limit, double upper_limit, int n) const; ///< Method to compute 1d numerical overlap.
    double gaussian(double x) const; ///< Method to compute Gaussian function value.
    double integrand(double x, const GaussianShell& shellA, const GaussianShell& shellB) const; ///< Method to compute the value of integrand.
    double computeSABx(const GaussianShell& shellB) const; ///< Method to compute SABx.
    double binomialCoefficient(int i) const; ///< Method to compute binomial coefficient.
    long long int factorial(int i) const; ///< Method to compute factorial.
    long long int doubleFactorial(int i) const; ///< Method to compute double factorial.
    double centerOfProductGaussian(const GaussianShell& shellB) const; ///< Method to compute center of product Gaussian.
};

#endif //PS_2_H
