#include "Convolution.h"

#include <algorithm>
#include <cmath>

namespace pals {

namespace {

constexpr double kInvSqrt2 = 0.70710678118654752440;
constexpr double kInvSqrt2Pi = 0.39894228040143267794;

// Phi(a) - Phi(b) for a > b, without cancellation in the tails
double normalProbability(double a, double b)
{
    if (b >= 0.0) return 0.5 * (std::erfc(b * kInvSqrt2) - std::erfc(a * kInvSqrt2));
    if (a <= 0.0) return 0.5 * (std::erfc(-a * kInvSqrt2) - std::erfc(-b * kInvSqrt2));
    return 1.0 - 0.5 * std::erfc(a * kInvSqrt2) - 0.5 * std::erfc(-b * kInvSqrt2);
}

double normalDensity(double z) { return kInvSqrt2Pi * std::exp(-0.5 * z * z); }

} // namespace

void GaussianConvolution::reset(double step)
{
    h_ = step;
    nodes_.clear();
}

double GaussianConvolution::operator()(double u, double sigma, double window) const
{
    if (nodes_.size() < 3) return 0.0;
    const double panel = 2.0 * h_;
    const long nPanels = static_cast<long>((nodes_.size() - 1) / 2);
    const double lo = u - window * sigma;
    const double hi = u + window * sigma;
    if (hi <= 0.0) return 0.0;
    const long m0 = std::max(0L, static_cast<long>(std::floor(lo / panel)));
    const long m1 = std::min(nPanels, static_cast<long>(std::ceil(hi / panel)));

    double sum = 0.0;
    // z of the left edge of the first panel; edges are shared between panels
    double a = (u - m0 * panel) / sigma;
    double phiA = normalDensity(a);
    for (long m = m0; m < m1; ++m) {
        const double xl = m * panel;
        const double b = (u - xl - panel) / sigma;     // z of the right edge (a > b)
        const double phiB = normalDensity(b);
        const double M0 = normalProbability(a, b);     // int_b^a phi
        const double M1 = phiB - phiA;                 // int_b^a z phi
        const double M2 = M0 + b * phiB - a * phiA;    // int_b^a z^2 phi

        const double dl = nodes_[2 * m], dm = nodes_[2 * m + 1], dr = nodes_[2 * m + 2];
        // D(x) = c0 + c1 (x - xm) + c2 (x - xm)^2,  x = u - sigma z,  v = u - xm
        const double c1 = (dr - dl) / panel;
        const double c2 = (dr - 2.0 * dm + dl) / (2.0 * h_ * h_);
        const double v = u - (xl + h_);
        const double A0 = dm + c1 * v + c2 * v * v;
        const double A1 = -sigma * (c1 + 2.0 * c2 * v);
        const double A2 = c2 * sigma * sigma;
        sum += A0 * M0 + A1 * M1 + A2 * M2;

        a = b;
        phiA = phiB;
    }
    return sum;
}

} // namespace pals
