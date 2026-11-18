#ifndef PALS_EMG_H
#define PALS_EMG_H

// Exponential decay (lifetime tau, starting at t0) convolved with a Gaussian
// (sigma), normalised to unit area:
//   f(t) = 1/(2 tau) exp(sigma^2/(2 tau^2) - (t - t0)/tau)
//          * erfc((sigma^2/tau - (t - t0)) / (sqrt(2) sigma))
// For z = (sigma/tau - (t-t0)/sigma)/sqrt(2) < 25 the direct form is finite
// (exp(...) <= exp(z^2)); further left exp(z^2) erfc(z) is taken from its
// asymptotic series (truncation error < 1e-16).
//
// The code is written once as a macro: it is compiled into the program and
// the same text is declared in the ROOT interpreter, so that it can be used
// in formulas from the config as pals_emg(x, sigma, t0, tau).
#define PALS_EMG_CODE                                                                          \
    inline double pals_emg(double t, double sigma, double t0, double tau)                    \
    {                                                                                          \
        const double u = (t - t0) / sigma;                                                     \
        const double a = sigma / tau;                                                          \
        const double z = (a - u) * 0.70710678118654752440;                                     \
        if (z < 25.0) return std::exp(a * (0.5 * a - u)) * std::erfc(z) / (2.0 * tau);         \
        const double h = 0.5 / (z * z);                                                        \
        const double s = 1.0 - h * (1.0 - 3.0 * h * (1.0 - 5.0 * h * (1.0 - 7.0 * h)));        \
        return std::exp(-0.5 * u * u) * s / (z * 1.77245385090551602730) / (2.0 * tau);        \
    }

#define PALS_STRINGIFY_(...) #__VA_ARGS__
#define PALS_STRINGIFY(...) PALS_STRINGIFY_(__VA_ARGS__)

#include <cmath>

PALS_EMG_CODE

// Source text of pals_emg for the ROOT interpreter.
inline const char* palsEmgSource() { return "#include <cmath>\n" PALS_STRINGIFY(PALS_EMG_CODE); }

#endif
