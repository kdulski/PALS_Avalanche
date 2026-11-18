#ifndef PALS_QUADRATURE_H
#define PALS_QUADRATURE_H

#include <cmath>
#include <vector>

namespace pals {

// Gauss-Hermite quadrature for the standard normal distribution:
//   E[f(Z)] = sum_k w_k f(x_k),  Z ~ N(0,1),  sum_k w_k = 1.
// Nodes from the Newton iteration on Hermite polynomials (Numerical Recipes gauher),
// converted from the weight exp(-x^2) to the normal density.
struct GaussHermite {
    std::vector<double> x, w;

    explicit GaussHermite(int n)
    {
        x.assign(n, 0.0);
        w.assign(n, 0.0);
        const double pim4 = 0.7511255444649425;   // pi^(-1/4)
        const int m = (n + 1) / 2;
        double z = 0.0, pp = 0.0;
        for (int i = 0; i < m; ++i) {
            if (i == 0) z = std::sqrt(2.0 * n + 1) - 1.85575 * std::pow(2.0 * n + 1, -0.16667);
            else if (i == 1) z -= 1.14 * std::pow(static_cast<double>(n), 0.426) / z;
            else if (i == 2) z = 1.86 * z - 0.86 * x[0];
            else if (i == 3) z = 1.91 * z - 0.91 * x[1];
            else z = 2.0 * z - x[i - 2];
            for (int it = 0; it < 100; ++it) {
                double p1 = pim4, p2 = 0.0;
                for (int j = 0; j < n; ++j) {
                    const double p3 = p2;
                    p2 = p1;
                    p1 = z * std::sqrt(2.0 / (j + 1)) * p2 - std::sqrt(static_cast<double>(j) / (j + 1)) * p3;
                }
                pp = std::sqrt(2.0 * n) * p2;
                const double z1 = z;
                z = z1 - p1 / pp;
                if (std::fabs(z - z1) <= 1e-14) break;
            }
            x[i] = z;
            x[n - 1 - i] = -z;
            w[i] = w[n - 1 - i] = 2.0 / (pp * pp);
        }
        for (int i = 0; i < n; ++i) {
            x[i] *= std::sqrt(2.0);
            w[i] /= std::sqrt(M_PI);
        }
    }
};

} // namespace pals

#endif
