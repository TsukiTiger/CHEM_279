#include "../inc/PS_2.h"
#include <iostream>
#include <cassert>

const double tolerance = 2e-1; // Tolerance for floating-point comparisons
const double integration_limit = 10.0; // Limit for numerical integration
const int integration_points = 10000; // Number of points for numerical integration

void runTest(const GaussianShell &shell1, const GaussianShell &shell2, const char *testName) {
    double analyticalResult = shell1.computeSABx(shell2);
    double numericalResult = shell1.computeNumericalOverlap_1d(shell2, -integration_limit, integration_limit, integration_points);

    std::cout << analyticalResult << std::endl;
    std::cout << numericalResult << std::endl;

    assert(fabs(analyticalResult - numericalResult) < tolerance);
    std::cout << testName << " Passed!" << std::endl;
}

int main() {

    // Part 1 Tests
    std::cout << "Running Part 1 Tests..." << std::endl;

    // 1. Two s-type functions both centered at the origin
    GaussianShell sShell1(0.0, 1.0, 0);
    runTest(sShell1, sShell1, "s-s centered at origin");

    // 2. An s-type function and a p-type function both centered at the origin
    GaussianShell pShell(0.0, 1.0, 1);
    runTest(sShell1, pShell, "s-p centered at origin");

    // 3. Repeat the above with an offset of 1
    GaussianShell offsetShell(1.0, 1.0, 0);
    runTest(offsetShell, sShell1, "s-s with an offset of 1");

    // Part 2 Tests
    std::cout << "Running Part 2 Tests..." << std::endl;

    // s-p shell pairs
    runTest(sShell1, pShell, "s-p shell pairs");

    // p-p shell pairs
    GaussianShell pShell2(0.0, 1.0, 1);
    runTest(pShell, pShell2, "p-p shell pairs");

    // Offsetting shells
    GaussianShell offsetPshell(1.0, 1.0, 1);
    runTest(offsetShell, offsetPshell, "s-p offsetting shells");

    // Check edge cases and other scenarios.
    // Example: Check the case where the angular momenta are mismatched, etc.

    // All tests passed
    std::cout << "All tests passed!" << std::endl;

    return 0;
}
