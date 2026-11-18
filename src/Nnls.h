#ifndef PALS_NNLS_H
#define PALS_NNLS_H

#include <cstddef>
#include <vector>

namespace pals {

// Dense column-major matrix (columns contiguous).
struct ColMatrix {
    std::size_t rows = 0, cols = 0;
    std::vector<double> a;
    ColMatrix() = default;
    ColMatrix(std::size_t r, std::size_t c) : rows(r), cols(c), a(r * c, 0.0) {}
    double& operator()(std::size_t i, std::size_t j) { return a[j * rows + i]; }
    double operator()(std::size_t i, std::size_t j) const { return a[j * rows + i]; }
    const double* col(std::size_t j) const { return a.data() + j * rows; }
};

// Non-negative least squares  min ||A x - b||, x >= 0
// Lawson & Hanson active-set algorithm ("Solving Least Squares Problems", 1974, ch. 23);
// the unconstrained subproblems on the passive set are solved by Householder QR.
struct NnlsResult {
    std::vector<double> x;
    double residualNorm = 0.0;
    int iterations = 0;
    bool converged = false;
};
NnlsResult nnls(const ColMatrix& A, const std::vector<double>& b, int maxIterations = 0);

// Solves the symmetric positive-definite system H x = g (H row-major n x n) by Cholesky;
// returns false if H is not positive definite.
bool choleskySolve(std::vector<double> H, std::vector<double>& x, const std::vector<double>& g, std::size_t n);

} // namespace pals

#endif
