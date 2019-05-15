# Parallel C++ linear algebra library with basic matrix operations 

In main.cpp file you could call any of the given functions you want to work with.
The set of functions is:
- matrix addition - matrix::add_parallel(M_1, M_2)
- matrix subtraction - matrix::sub_parallel(M_1, M_2);
- matrix multiplication - matrix::mul_parallel(M_1, M_2);
- inverse of the matrix - matrix::inverse(M).

The Instruction:
1. At first you have to initialize the matrix(e.g. matrix M_1(3, 3);) or matrices. 
2. Assign values to your initialized matrix or matrices.
3. Call one of the functions in the above set.
4. Use M.to_string() to have a nice output of your M matrix.

A couple of the examples are situated right inside the main.cpp, so as you could analogically make other calculations. 
For the source code look into matrix.cpp and matrix.hpp files.

© Created by Shtokhman Yuriy, Taras Hrytsiv and Anna Hayda.
