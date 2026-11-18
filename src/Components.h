#ifndef PALS_COMPONENTS_H
#define PALS_COMPONENTS_H

#include <string>

namespace pals {

// How a single fit parameter behaves during the iterative fit.
//   Free    - free within its global limits
//   Fixed   - never released
//   Partial - released within +- (partial variation) of its current value;
//             the window is re-centred in every iteration
enum class ParMode { Free, Fixed, Partial };

inline const char* toString(ParMode m)
{
    switch (m) {
        case ParMode::Free:    return "free";
        case ParMode::Fixed:   return "fixed";
        case ParMode::Partial: return "partial";
    }
    return "?";
}

// Gaussian component of the instrumental resolution function.
struct ResolutionComponent {
    double sigma = 0.0;          // [ns]
    bool sigmaAuto = false;      // estimate sigma from the left half-width of the peak
    double intensity = 1.0;      // fraction of the resolution function
    ParMode sigmaMode = ParMode::Free;
    ParMode intensityMode = ParMode::Free;
};

// Exponential (lifetime) component.
// Intensities with intensityMode == Free form a group whose relative
// fractions are described by the hyperspherical parametrisation and which
// shares whatever is left after subtracting the Fixed/Partial ones.
// For the p-Ps component (pPs == true) "intensity" is the ratio
// I(p-Ps)/I(o-Ps) (1/3 in vacuum); it is always kept fixed.
struct LifetimeComponent {
    double tau = 0.0;            // [ns]
    double intensity = 0.0;
    ParMode tauMode = ParMode::Free;
    ParMode intensityMode = ParMode::Free;
    bool pPs = false;
    // optional log-normal distribution of tau around the mean "tau":
    // width = standard deviation of ln(tau) (0 = one exponential)
    bool hasWidth = false;
    double width = 0.0;
    ParMode widthMode = ParMode::Free;
};

} // namespace pals

#endif
