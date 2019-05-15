#include <iostream>
#include "matrix.hpp"


int main() {
    // Some examples of usage
    // For more detailed instructions open README.md file
    matrix M_1(3, 3);
    std::vector<std::vector<double>> values_1 = {{1, 2, 3},
                                                 {4, 5, 9},
                                                 {3, 3, 3}};
    M_1 = values_1;

    matrix M_2(3, 3);
    std::vector<std::vector<double>> values_2 = {{4, 2, 1},
                                                 {4, 8, 3},
                                                 {2, 2, 2}};
    M_2 = values_2;

    matrix M_inv(3, 3);
    std::vector<std::vector<double>> values_inv = {{1, 2, 3},
                                                   {4, 5, 9},
                                                   {3, 3, 3}};
    M_inv = values_inv;


    matrix M_add_parallel(3, 3);
    M_add_parallel = matrix::add_parallel(M_1, M_2);
    std::cout << M_add_parallel.to_string();

    matrix M_inverse(3, 3);
    M_inverse = matrix::inverse(M_inv);
    std::cout << M_inverse.to_string() << "\n";

    return 0;

}