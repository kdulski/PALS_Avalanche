#ifndef PALS_INTENSITIES_H
#define PALS_INTENSITIES_H

// Hyperspherical parametrisation of n fractions w_1..w_n with sum(w) = 1.
//
// One parameter per component (the same layout as in the original code):
//   q[0]            = r                 in [0, 1]
//   q[1] .. q[n-2]  = angles phi_k      in [0, pi/2]
//   q[n-1]          = not used (must be kept fixed in the fit)
//
//   x_1 = cos(phi_1)
//   x_k = sin(phi_1)...sin(phi_{k-1}) cos(phi_k)       k = 2..n-2
//   x_{n-1} = sin(phi_1)...sin(phi_{n-2})
//   w_k = r * x_k / sum(x)     k = 1..n-1
//   w_n = 1 - r
//
// Every point of the box r in [0,1], phi in [0,pi/2] gives w_k >= 0 with
// sum(w) = 1, so no constraint has to be imposed on Minuit. This is the
// generic form of the explicit formulas written case-by-case (n <= 10) in
// GetIntensityParameterNew / SetIntensityParameter.

#include <cmath>
#include <cstddef>
#include <numeric>
#include <vector>

namespace pals {

// q -> w ; q and w have n elements
inline void hypersphericalToFractions(const double* q, std::size_t n, double* w)
{
    if (n == 0) return;
    if (n == 1) { w[0] = 1.0; return; }
    const std::size_t m = n - 1;          // components sharing r
    const double r = q[0];
    double sinProduct = 1.0;
    double sum = 0.0;
    for (std::size_t k = 0; k + 1 < m; ++k) {
        w[k] = sinProduct * std::cos(q[1 + k]);
        sinProduct *= std::sin(q[1 + k]);
        sum += w[k];
    }
    w[m - 1] = sinProduct;
    sum += sinProduct;
    const double scale = sum > 0.0 ? r / sum : 0.0;
    for (std::size_t k = 0; k < m; ++k) w[k] *= scale;
    w[m] = 1.0 - r;
}

// w -> q (used for the initial values); w does not have to be normalised.
// The unused last parameter is set to 0.
inline void fractionsToHyperspherical(const double* wIn, std::size_t n, double* q)
{
    if (n == 0) return;
    for (std::size_t k = 0; k < n; ++k) q[k] = 0.0;
    if (n == 1) return;

    double total = 0.0;
    for (std::size_t k = 0; k < n; ++k) total += wIn[k] > 0.0 ? wIn[k] : 0.0;
    std::vector<double> w(n, 1.0 / n);
    if (total > 0.0)
        for (std::size_t k = 0; k < n; ++k) w[k] = (wIn[k] > 0.0 ? wIn[k] : 0.0) / total;

    const std::size_t m = n - 1;
    q[0] = 1.0 - w[m];
    if (m < 2) return;

    // x_k proportional to w_k (normalisation of x does not matter for the angles)
    // phi_k = atan2( |(x_{k+1}, ..., x_{m})|, x_k )
    std::vector<double> tailNorm(m + 1, 0.0);
    for (std::size_t k = m; k-- > 0;)
        tailNorm[k] = std::hypot(tailNorm[k + 1], w[k]);
    for (std::size_t k = 0; k + 1 < m; ++k) {
        const double tail = tailNorm[k + 1];
        q[1 + k] = (tail == 0.0 && w[k] == 0.0) ? 0.25 * M_PI : std::atan2(tail, w[k]);
    }
}

// Fractions of a group of components: members with isFree == false take their
// parameter as an absolute fraction, the free ones share the rest (1 - sum of
// absolute ones, at least 0) according to the hyperspherical parametrisation
// of their parameters (in member order).
inline std::vector<double> groupFractions(const std::vector<double>& pars, const std::vector<bool>& isFree)
{
    std::vector<double> out(pars.size(), 0.0);
    std::vector<double> q;
    std::vector<std::size_t> freeIdx;
    double absSum = 0.0;
    for (std::size_t k = 0; k < pars.size(); ++k) {
        if (isFree[k]) { freeIdx.push_back(k); q.push_back(pars[k]); }
        else { out[k] = pars[k]; absSum += pars[k]; }
    }
    if (freeIdx.empty()) return out;
    std::vector<double> w(q.size());
    hypersphericalToFractions(q.data(), q.size(), w.data());
    const double rest = absSum < 1.0 ? 1.0 - absSum : 0.0;
    for (std::size_t k = 0; k < freeIdx.size(); ++k) out[freeIdx[k]] = rest * w[k];
    return out;
}

} // namespace pals

#endif
