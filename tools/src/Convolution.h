#ifndef PALS_CONVOLUTION_H
#define PALS_CONVOLUTION_H

#include <cstddef>
#include <vector>

namespace pals {

// Numerical convolution of a function D(x), defined for x >= 0 (zero below),
// with a unit-area Gaussian:
//   C(u) = integral_0^X D(x) G(u - x; sigma) dx
//
// D is given on a uniform grid x_k = k h (an odd number of nodes, so that the
// grid consists of panels [x_2m, x_2m+2]). On every panel D is replaced by the
// quadratic through its three nodes and the product with the Gaussian is
// integrated exactly (moments of the normal distribution), so
//   - the step at x = 0 is treated exactly,
//   - the error is O(h^3) per panel and does not depend on sigma/h aliasing.
// Only panels within +- window*sigma of u contribute.
class GaussianConvolution {
public:
    void reset(double step);
    void append(double value) { nodes_.push_back(value); }

    double step() const { return h_; }
    std::size_t size() const { return nodes_.size(); }
    // end of the grid (valid when the number of nodes is odd)
    double range() const { return nodes_.empty() ? 0.0 : (nodes_.size() - 1) * h_; }

    double operator()(double u, double sigma, double window) const;

private:
    double h_ = 0.0;
    std::vector<double> nodes_;
};

} // namespace pals

#endif
