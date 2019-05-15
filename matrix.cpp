#include "matrix.hpp"
#include <thread>
#include <mutex>
#include <iostream>


std::mutex f_lock;


std::vector<size_t> make_vector_steps_for_threads(size_t num, size_t num_threads){
    int step = static_cast<int>(num/num_threads);
    std::vector <size_t> steps;
    for (size_t i = 0; i < num_threads; i++) steps.push_back(i*step);
    steps.push_back(num);
    return steps;
}

void subtract_matrix(matrix &addition_matrix, const matrix &m1, const matrix &m2, size_t from, size_t to) {
    for (size_t i = from; i < to; i++) {
        for (size_t j = 0; j < m1.get_cols(); j++) {
            f_lock.lock();
            addition_matrix(i, j) = m1(i, j) - m2(i, j);
            f_lock.unlock();
        }
    }
}

void add_matrix(matrix &addition_matrix, const matrix &m1, const matrix &m2, size_t from, size_t to) {
    for (size_t i = from; i < to; i++) {
        for (size_t j = 0; j < m1.get_cols(); j++) {
            f_lock.lock();
            addition_matrix(i, j) = m1(i, j) + m2(i, j);
            f_lock.unlock();
        }
    }
}

void mul_matrix(matrix &mull_matrix, const matrix &m1, const matrix &m2, size_t from, size_t to){
    for (size_t i = from; i < to; i++) {
        for (size_t j = 0; j < m2.get_cols(); j++) {
            for (size_t k = 0; k < m2.get_rows(); k++) {
                f_lock.lock();
                mull_matrix(i, j) += m1(i, k) * m2(k, j);
                f_lock.unlock();
            }
        }
    }
}

matrix matrix::add(const matrix &m1, const matrix &m2) {
    if (m1.get_cols() != m2.get_cols() || m1.get_rows() != m2.get_rows()) {
        std::__throw_invalid_argument("matrices are of different sizes");
    }
    matrix addition_matrix(m1.get_rows(), m1.get_cols());
    for (size_t i = 0; i < m1.get_rows(); i++) {
        for (size_t j = 0; j < m1.get_cols(); j++) {
            addition_matrix(i, j) = m1(i, j) + m2(i, j);
        }
    }
    return addition_matrix;
}

matrix matrix::sub(const matrix &m1, const matrix &m2) {
    // check matrix sizes
    if (m1.get_cols() != m2.get_cols() || m1.get_rows() != m2.get_rows()) {
        std::__throw_invalid_argument("matrices are of different sizes");
    }
    matrix addition_matrix(m1.get_rows(), m1.get_cols());
    for (size_t i = 0; i < m1.get_rows(); i++) {
        for (size_t j = 0; j < m1.get_cols(); j++) {
            addition_matrix(i, j) = m1(i, j) - m2(i, j);
        }
    }
    return addition_matrix;
}

matrix matrix::mul(const matrix &m1,const matrix &m2) {
    if (m1.get_cols() != m2.get_rows() || m1.get_rows() != m2.get_cols()) {
        std::__throw_invalid_argument("matrices are of different sizes");
    }

    matrix multiplication_matrix(m1.get_rows(), m1.get_rows());

    for (size_t i = 0; i < m1.get_rows(); i++) {
        for (size_t j = 0; j < m2.get_cols(); j++) {
            for (size_t k = 0; k < m2.get_rows(); k++) {
                multiplication_matrix(i, j) += m1(i, k) * m2(k, j);
            }
        }
    }

    return multiplication_matrix;
}

void paralleling(const matrix &m1, matrix &U, matrix &L, matrix &identity, size_t from, size_t to){
    f_lock.lock();
    size_t dim = m1.get_rows();
    f_lock.unlock();

    for (size_t i = from; i < to; ++i) {
        for (size_t k = i; k < dim; ++k) {
            double sum = 0;
            for (size_t j = 0; j < i; ++j) {
                f_lock.lock();
                sum += (L(i, j) * U(j, k));
                f_lock.unlock();
            }
            f_lock.lock();
            U(i, k) = m1(i, k) - sum;
            f_lock.unlock();
        }
        for (size_t k = i; k < dim; ++k) {
            if (i == k) {
                f_lock.lock();
                L(i, i) = 1;
                f_lock.unlock();
            } else {
                double sum = 0;
                for (size_t j = 0; j < i; ++j) {
                    f_lock.lock();
                    sum += (L(k, j) * U(j, i));
                    f_lock.unlock();
                }
                f_lock.lock();
                L(k, i) = (m1(k, i) - sum) / U(i, i);
                f_lock.unlock();
            }
        }
    }

    for (size_t i = from; i < to; ++i) {
        for (size_t j = 0; j < dim; ++j) {
            if (i == j) {
                f_lock.lock();
                identity(i, j) = 1;
                f_lock.unlock();
            }
        }
    }
}

void inverse_m(const matrix &m1, matrix &U, matrix &L, matrix &identity,size_t from, size_t to ) {
    f_lock.lock();
    size_t dim = m1.get_rows();
    f_lock.unlock();

    std::vector<double> z(dim, 0);
    std::vector<double> b(dim, 0);
    paralleling(m1, U, L, identity, from, to);
}

matrix matrix::inverse_parallel(const matrix &m1) {
    if (m1.get_cols() != m1.get_rows()) {
        std::__throw_invalid_argument("matrix is not square");
    }

    size_t num_threads = std::thread::hardware_concurrency();
    if (num_threads > m1.get_rows()) num_threads = m1.get_rows();

    size_t dim = m1.get_rows();
    matrix L(dim, dim);
    matrix U(dim, dim);
    matrix prefinal(dim, dim);
    matrix inverted(dim, dim);

    matrix identity(dim, dim);

    std::vector<size_t> steps = make_vector_steps_for_threads(m1.get_rows(), num_threads);

    std::thread myThreads[num_threads];
    for (int i = 0; i < num_threads; i++) {
        myThreads[i] = std::thread(inverse_m, std::ref(m1), std::ref(U), std::ref(L),  std::ref(identity), steps[i], steps[i+1]);
    }
    for (int i = 0; i < num_threads; i++) myThreads[i].join();

    double sum;
    double z[dim];
    double b[dim];
    for (size_t c = 0; c < dim; ++c) {
        for (size_t i = 0; i < dim; ++i) {

            sum = identity(c, i);
            for (size_t j = 0; j < i; ++j) {
                sum -= L(i, j) * z[j];
            }
            z[i] = sum / L(i, i);
        }
        for (int i = dim-1; i >= 0; --i) {
            sum = z[i];
            for (size_t j = i+1; j < dim; ++j) {
                sum = sum - U(i, j)*b[j];
            }
            b[i] = sum / U(i, i);
            prefinal(c, i) = b[i];

        }
    }

    for (size_t j = 0; j < dim; ++j) {
        for (size_t i = 0; i < dim; ++i) {
            inverted(j, i) = prefinal(i, j);
        }
    }

    return inverted;
}

matrix matrix::inverse(const matrix &m1) {
    if (m1.get_cols() != m1.get_rows()) {
        std::__throw_invalid_argument("matrix is not square");
    }
    size_t dim = m1.get_rows();
    std::cout << dim;
    matrix L(dim, dim);
    matrix U(dim, dim);
    matrix identity(dim, dim);
    matrix prefinal(dim, dim);
    matrix inverted(dim, dim);
    double sum;

    for (size_t i = 0; i < dim; ++i) {
        for (size_t k = i; k < dim; ++k) {
            sum = 0;
            for (size_t j = 0; j < i; ++j) {
                sum += (L(i, j) * U(j, k));
            }
            U(i, k) = m1(i, k) - sum;
        }
        for (size_t k = i; k < dim; ++k) {
            if (i == k) {
                L(i, i) = 1;
            } else {
                sum = 0;
                for (size_t j = 0; j < i; ++j) {
                    sum += (L(k, j) * U(j, i));
                }
                L(k, i) = (m1(k, i) - sum) / U(i, i);
            }
        }
    }

    for (size_t i = 0; i < dim; ++i) {
        for (size_t j = 0; j < dim; ++j) {
            if (i == j) { identity(i, j) = 1; }
        }
    }

    double z[dim];
    double b[dim];
    for (size_t c = 0; c < dim; ++c) {
        for (size_t i = 0; i < dim; ++i) {

            sum = identity(c, i);
            for (size_t j = 0; j < i; ++j) {
                sum -= L(i, j) * z[j];
            }
        z[i] = sum / L(i, i);
        }
        for (int i = dim-1; i >= 0; --i) {
            sum = z[i];
            for (size_t j = i+1; j < dim; ++j) {
                sum = sum - U(i, j)*b[j];
            }
            b[i] = sum / U(i, i);
            prefinal(c, i) = b[i];

        }
    }

    for (size_t j = 0; j < dim; ++j) {
        for (size_t i = 0; i < dim; ++i) {
            inverted(j, i) = prefinal(i, j);
        }
    }

    return inverted;
}

matrix matrix::sub_parallel(const matrix &m1, const matrix &m2) {
    size_t num_threads = std::thread::hardware_concurrency();
    if (num_threads > m1.get_rows()) num_threads = m1.get_rows();
    if (m1.get_cols() != m2.get_cols() || m1.get_rows() != m2.get_rows()) {
        std::__throw_invalid_argument("matrices are of different sizes");
    }
    matrix addition_matrix(m1.get_rows(), m1.get_cols());
    std::vector<size_t> steps = make_vector_steps_for_threads(m1.get_rows(), num_threads);
    std::thread myThreads[num_threads];
    for (int i = 0; i < num_threads; i++) {
        myThreads[i] = std::thread(subtract_matrix, std::ref(addition_matrix), std::ref(m1), std::ref(m2), steps[i], steps[i + 1]);
    }
    for (int i = 0; i < num_threads; i++) myThreads[i].join();

    return addition_matrix;

}

matrix matrix::add_parallel(const matrix &m1, const matrix &m2) {
    size_t num_threads = std::thread::hardware_concurrency();
    if (num_threads > m1.get_rows()) num_threads = m1.get_rows();

    if (m1.get_cols() != m2.get_cols() || m1.get_rows() != m2.get_rows()) {
        std::__throw_invalid_argument("matrices are of different sizes");
    }
    matrix addition_matrix(m1.get_rows(), m1.get_cols());
    std::vector<size_t> steps = make_vector_steps_for_threads(m1.get_rows(), num_threads);
    std::thread myThreads[num_threads];
    for (int i = 0; i < num_threads; i++) {
        myThreads[i] = std::thread(add_matrix, std::ref(addition_matrix), std::ref(m1), std::ref(m2), steps[i], steps[i + 1]);
    }
    for (int i = 0; i < num_threads; i++) myThreads[i].join();

    return addition_matrix;

}

matrix matrix::mul_parallel(const matrix &m1, const matrix &m2) {
    size_t num_threads = std::thread::hardware_concurrency();
    if (num_threads > m1.get_rows()) num_threads = m1.get_rows();

    if (m1.get_cols() != m2.get_rows() || m1.get_rows() != m2.get_cols()) {
        std::__throw_invalid_argument("matrices are of different sizes");
    }

    matrix multiplication_matrix(m1.get_rows(), m1.get_rows());
    std::vector<size_t> steps = make_vector_steps_for_threads(m1.get_rows(), num_threads);
    std::thread myThreads[num_threads];

    for (int i = 0; i < num_threads; i++) {
        myThreads[i] = std::thread(mul_matrix, std::ref(multiplication_matrix), std::ref(m1), std::ref(m2), steps[i],
                                   steps[i + 1]);
    }

    for (int i = 0; i < num_threads; i++) myThreads[i].join();

    return multiplication_matrix;
}

matrix &matrix::operator=(const std::vector<std::vector<double>> &init_values) {
    if (init_values.size() != rows || init_values[0].size() != cols) {
        std::__throw_invalid_argument("size of init 2d array doesn't match matrix dimensions");
    }
    for (size_t i = 0; i < init_values.size(); ++i) {
        for (size_t j = 0; j < init_values[i].size(); ++j) {
            data[i][j] = init_values[i][j];
        }
    }
    return *this;
}
