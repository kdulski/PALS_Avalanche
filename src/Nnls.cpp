#include "Nnls.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace pals {

namespace {

// least squares on the columns 'cols' of A by Householder QR; returns z (size = cols.size())
std::vector<double> lstsq(const ColMatrix& A, const std::vector<double>& b, const std::vector<std::size_t>& cols)
{
    const std::size_t m = A.rows, p = cols.size();
    std::vector<double> R(m * p);                 // column-major copy
    for (std::size_t k = 0; k < p; ++k) std::copy(A.col(cols[k]), A.col(cols[k]) + m, R.begin() + k * m);
    std::vector<double> y = b;
    std::vector<double> diag(p, 0.0);
    for (std::size_t k = 0; k < p && k < m; ++k) {
        double* c = R.data() + k * m;
        double norm = 0.0;
        for (std::size_t i = k; i < m; ++i) norm += c[i] * c[i];
        norm = std::sqrt(norm);
        if (norm == 0.0) { diag[k] = 0.0; continue; }
        const double alpha = c[k] > 0 ? -norm : norm;
        // v = c[k..] - alpha e_k, stored in c[k..]
        c[k] -= alpha;
        double vnorm2 = 0.0;
        for (std::size_t i = k; i < m; ++i) vnorm2 += c[i] * c[i];
        diag[k] = alpha;
        if (vnorm2 == 0.0) continue;
        for (std::size_t j = k + 1; j < p; ++j) {
            double* d = R.data() + j * m;
            double s = 0.0;
            for (std::size_t i = k; i < m; ++i) s += c[i] * d[i];
            s = 2.0 * s / vnorm2;
            for (std::size_t i = k; i < m; ++i) d[i] -= s * c[i];
        }
        double s = 0.0;
        for (std::size_t i = k; i < m; ++i) s += c[i] * y[i];
        s = 2.0 * s / vnorm2;
        for (std::size_t i = k; i < m; ++i) y[i] -= s * c[i];
    }
    // back substitution R z = (Q^T b)[0..p)
    std::vector<double> z(p, 0.0);
    double rmax = 0.0;
    for (std::size_t k = 0; k < p; ++k) rmax = std::max(rmax, std::fabs(diag[k]));
    for (std::size_t k = p; k-- > 0;) {
        if (std::fabs(diag[k]) <= 1e-13 * rmax) { z[k] = 0.0; continue; }   // rank deficient: drop
        double s = y[k];
        for (std::size_t j = k + 1; j < p; ++j) s -= R[j * m + k] * z[j];
        z[k] = s / diag[k];
    }
    return z;
}

} // namespace

NnlsResult nnls(const ColMatrix& A, const std::vector<double>& b, int maxIterations)
{
    const std::size_t m = A.rows, n = A.cols;
    if (maxIterations <= 0) maxIterations = static_cast<int>(3 * n);
    NnlsResult res;
    res.x.assign(n, 0.0);
    std::vector<bool> passive(n, false);
    std::vector<double> r = b, w(n);

    double anorm = 0.0;
    for (std::size_t j = 0; j < n; ++j) {
        double s = 0.0;
        for (std::size_t i = 0; i < m; ++i) s += std::fabs(A(i, j));
        anorm = std::max(anorm, s);
    }
    const double tol = 10.0 * std::numeric_limits<double>::epsilon() * anorm * std::max(m, n);

    auto gradient = [&]() {
        for (std::size_t j = 0; j < n; ++j) {
            const double* c = A.col(j);
            double s = 0.0;
            for (std::size_t i = 0; i < m; ++i) s += c[i] * r[i];
            w[j] = s;
        }
    };
    auto residual = [&]() {
        r = b;
        for (std::size_t j = 0; j < n; ++j) {
            if (res.x[j] == 0.0) continue;
            const double* c = A.col(j);
            for (std::size_t i = 0; i < m; ++i) r[i] -= c[i] * res.x[j];
        }
    };

    gradient();
    int it = 0;
    while (it < maxIterations) {
        // most positive gradient among the active (zero) variables
        std::size_t t = n;
        double wmax = tol;
        for (std::size_t j = 0; j < n; ++j)
            if (!passive[j] && w[j] > wmax) { wmax = w[j]; t = j; }
        if (t == n) { res.converged = true; break; }
        passive[t] = true;

        // inner loop: keep the passive solution feasible
        while (true) {
            ++it;
            std::vector<std::size_t> P;
            for (std::size_t j = 0; j < n; ++j) if (passive[j]) P.push_back(j);
            const std::vector<double> z = lstsq(A, b, P);
            bool feasible = true;
            for (std::size_t k = 0; k < P.size(); ++k) if (z[k] <= 0.0) { feasible = false; break; }
            if (feasible) {
                for (std::size_t k = 0; k < P.size(); ++k) res.x[P[k]] = z[k];
                break;
            }
            double alpha = 1.0;
            for (std::size_t k = 0; k < P.size(); ++k) {
                if (z[k] > 0.0) continue;
                const double xj = res.x[P[k]];
                const double a = xj / (xj - z[k]);
                alpha = std::min(alpha, a);
            }
            for (std::size_t k = 0; k < P.size(); ++k) {
                double& xj = res.x[P[k]];
                xj += alpha * (z[k] - xj);
                if (xj <= tol) { xj = 0.0; passive[P[k]] = false; }
            }
            if (it >= maxIterations) break;
        }
        residual();
        gradient();
    }
    residual();
    double rn = 0.0;
    for (double v : r) rn += v * v;
    res.residualNorm = std::sqrt(rn);
    res.iterations = it;
    return res;
}

bool choleskySolve(std::vector<double> H, std::vector<double>& x, const std::vector<double>& g, std::size_t n)
{
    for (std::size_t j = 0; j < n; ++j) {
        double d = H[j * n + j];
        for (std::size_t k = 0; k < j; ++k) d -= H[j * n + k] * H[j * n + k];
        if (d <= 0.0) return false;
        d = std::sqrt(d);
        H[j * n + j] = d;
        for (std::size_t i = j + 1; i < n; ++i) {
            double s = H[i * n + j];
            for (std::size_t k = 0; k < j; ++k) s -= H[i * n + k] * H[j * n + k];
            H[i * n + j] = s / d;
        }
    }
    x = g;
    for (std::size_t i = 0; i < n; ++i) {
        double s = x[i];
        for (std::size_t k = 0; k < i; ++k) s -= H[i * n + k] * x[k];
        x[i] = s / H[i * n + i];
    }
    for (std::size_t i = n; i-- > 0;) {
        double s = x[i];
        for (std::size_t k = i + 1; k < n; ++k) s -= H[k * n + i] * x[k];
        x[i] = s / H[i * n + i];
    }
    return true;
}

} // namespace pals
