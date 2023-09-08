#include <armadillo>
#include <iostream>

int main() {
    arma::vec2 test_vec;  // Corrected to vec2
    test_vec(0) = 3;
    test_vec(1) = 6;
    std::cout << test_vec << std::endl;

    return 0;
}
