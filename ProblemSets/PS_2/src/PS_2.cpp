#include "../inc/PS_2.h"

GaussianShell::GaussianShell() : center(0.0), exponent(0.0), angularMomentum(0) {}

GaussianShell::GaussianShell(double c, double e, int l) : center(c), exponent(e), angularMomentum(l) {}

double GaussianShell::computeAnalyticalOverlap_1d(const GaussianShell& shellB) const {
    double alphaPlusBeta = exponent + shellB.exponent;
    double productCenter = centerOfProductGaussian(shellB);
    double prefactor = exp(-exponent * shellB.exponent * pow(center - shellB.center, 2) / alphaPlusBeta) * sqrt(M_PI / alphaPlusBeta);

    double sum = 0.0;
    for (int i = 0; i <= angularMomentum; ++i) {
        for (int j = 0; j <= shellB.angularMomentum; ++j) {
            sum += binomialCoefficient(i) * shellB.binomialCoefficient(j) * doubleFactorial(i + j - 1)
                   * pow(productCenter - center, angularMomentum - i)
                   * pow(productCenter - shellB.center, shellB.angularMomentum - j)
                   / pow(2.0 * alphaPlusBeta, (i + j) / 2.0);
        }
    }
    return prefactor * sum;
}

double GaussianShell::computeNumericalOverlap_1d(const GaussianShell& shellB, double lower_limit, double upper_limit, int n) const {
    double h = (upper_limit - lower_limit) / n;
    double sum = 0.0;

    // Evaluate the function at the endpoints
    sum += integrand(lower_limit, *this, shellB);
    sum += integrand(upper_limit, *this, shellB);

    // Evaluate the function at the interior points
    for (int i = 1; i < n; ++i) {
        double x = lower_limit + i * h;
        double f = integrand(x, *this, shellB);
        sum += 2 * f;
    }

    return sum * h / 2.0;
}

double GaussianShell::computeSABx(const GaussianShell& shellB) const {
    double alpha = this->exponent;
    double beta = shellB.exponent;
    double XA = this->center;
    double XB = shellB.center;
    int lA = this->angularMomentum;
    int lB = shellB.angularMomentum;

    double XP = centerOfProductGaussian(shellB);
    double coeff = exp(-alpha * beta * pow(XA - XB, 2) / (alpha + beta)) * sqrt(M_PI / (alpha + beta));

    double sum = 0.0;

    for(int i = 0; i <= lA; i++) {
        for(int j = 0; j <= lB; j++) {
            double binomCoeffA = this->binomialCoefficient(i);
            double binomCoeffB = shellB.binomialCoefficient(j);
            double term1 = binomCoeffA * binomCoeffB * pow(XP - XA, lA - i) * pow(XP - XB, lB - j);
            double term2 = pow(2 * (alpha + beta), (i + j) / 2.0);
            double term3 = (i + j - 1 >= 0) ? this->doubleFactorial(i + j - 1) : 1;

            sum += term1 * term3 / term2;
        }
    }

    return coeff * sum;
}

double GaussianShell::gaussian(double x) const {
    return pow(x - center, angularMomentum) * exp(-exponent * pow(x - center, 2));
}

double GaussianShell::integrand(double x, const GaussianShell& shellA, const GaussianShell& shellB) const{
    return fabs(shellA.gaussian(x) * shellB.gaussian(x));
}

double GaussianShell::binomialCoefficient(int i) const {
    if (i < 0 || i > angularMomentum) return 0.0; // Return 0 for invalid i
    return static_cast<double>(factorial(angularMomentum)) /
           (static_cast<double>(factorial(i)) * static_cast<double>(factorial(angularMomentum - i)));
}

long long int GaussianShell::factorial(int i) const{
    if (i < 0) return 0; // Return 0 for negative numbers as they don't have a factorial
    long long int result = 1;
    for (int j = 2; j <= i; j++) {
        result *= j;
    }
    return result;
}

long long int GaussianShell::doubleFactorial(int i) const{
    if (i <= 0) return 1;
    long long int result = 1;
    for (int j = i; j > 0; j -= 2) {
        result *= j;
    }
    return result;
}

double GaussianShell::centerOfProductGaussian(const GaussianShell& shellB) const {
    return (exponent * center + shellB.exponent * shellB.center) / (exponent + shellB.exponent);
}
