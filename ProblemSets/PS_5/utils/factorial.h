#ifndef UTILS_FACTORIAL_H
#define UTILS_FACTORIAL_H

/**
 * @file factorial.h
 * @brief Contains utility functions related to factorials and binomial coefficients.
 */

/**
 * Computes the factorial of a given integer.
 *
 * @param n Integer whose factorial needs to be computed.
 * @return Returns n! (n factorial).
 */
inline double factorial(int n) {
    int result = 1;
    for (int i = n; i > 0; i--) {
        result *= i;
    }
    return result;
}

/**
 * Computes the binomial coefficient for a pair of integers.
 *
 * @note The binomial coefficient is defined only when m >= n.
 *
 * @param m The top number in the binomial coefficient.
 * @param n The bottom number in the binomial coefficient.
 * @return Returns the binomial coefficient C(m, n).
 */
inline double binomialCoef(int m, int n) {
    if (m < n) {
        return 0;
    }
    return factorial(m) / (factorial(n) * factorial(m - n));
}

/**
 * Computes the double factorial of a given integer.
 *
 * @note The double factorial of n is the product of all integers from 1 to n that have the same parity (odd/even) as n.
 *
 * @param n Integer whose double factorial needs to be computed.
 * @return Returns the double factorial of n.
 */
inline double doubleFactorial(int n) {
    int result = 1;
    for (int i = n; i > 0; i -= 2) {
        result *= i;
    }
    return result;
}

#endif // UTILS_FACTORIAL_H
