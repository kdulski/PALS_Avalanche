#ifndef PALS_LIFETIMEMODEL_H
#define PALS_LIFETIMEMODEL_H

#include "Components.h"
#include "Resolution.h"

#include <string>
#include <vector>

namespace pals {

// Exponential decay convolved with a Gaussian, unit area (see Emg.h).
double emg(double t, double sigma, double t0, double tau);

// Model of a PAL spectrum:
//   N(t) = B + A * binWidth * sum_i w_i sum_j I_j f(t; sigma_i, t0_i, tau_j)
// A is the number of annihilation events (counts), so its initial value is
// simply the background-subtracted number of counts.
//
// Parameter layout:
//   0                    background B   [counts / bin]
//   1                    area A         [counts]
//   2 + 3i               sigma_i        (resolution component i)
//   3 + 3i               t0_i           (offset of resolution component i)
//   4 + 3i               fraction parameter of resolution component i
//   2 + 3R + 2j          tau_j          (lifetime component j; mean lifetime)
//   3 + 3R + 2j          intensity parameter of lifetime component j
//   2 + 3R + 2L + k      width (sd of ln tau) of the k-th component with a width
//
// A component with a width is a log-normal distribution of lifetimes with mean
// tau_j, integrated with Gauss-Hermite quadrature.
//
// Intensity parameter meaning:
//   Free     -> hyperspherical parameter of the "free" group (see Intensities.h)
//   Fixed /
//   Partial  -> absolute intensity
//   p-Ps     -> ratio I(p-Ps)/I(o-Ps) (fixed)
//
// The model can be passed directly to TF1 (functor with operator()).
class LifetimeModel {
public:
    LifetimeModel(std::vector<ResolutionComponent> resolution,
                  std::vector<LifetimeComponent> lifetimes,
                  double binWidth, double oPsThreshold);

    static constexpr int kBackground = 0;
    static constexpr int kArea = 1;

    int nResolution() const { return static_cast<int>(resolution_.size()); }
    int nLifetimes() const { return static_cast<int>(lifetimes_.size()); }
    int nPar() const { return 2 + 3 * nResolution() + 2 * nLifetimes() + nWidths_; }

    int sigmaPar(int i) const { return 2 + 3 * i; }
    int offsetPar(int i) const { return 3 + 3 * i; }
    int resolutionFractionPar(int i) const { return 4 + 3 * i; }
    int tauPar(int j) const { return 2 + 3 * nResolution() + 2 * j; }
    int intensityPar(int j) const { return 3 + 3 * nResolution() + 2 * j; }
    // parameter of the width of component j, -1 if it has none
    int widthPar(int j) const { return widthPar_[j]; }

    const std::vector<ResolutionComponent>& resolution() const { return resolution_; }
    const std::vector<LifetimeComponent>& lifetimes() const { return lifetimes_; }
    int pPsIndex() const { return pPs_; }
    double oPsThreshold() const { return oPsThreshold_; }
    double binWidth() const { return binWidth_; }

    // Parameter indices of the hyperspherical groups, in group order:
    // first = r, then angles, last = unused (to be fixed). A group of one
    // component consists only of the unused parameter.
    std::vector<int> resolutionGroupPars() const { return resLayout_.groupPars(); }
    const ResolutionLayout& resolutionLayout() const { return resLayout_; }
    std::vector<int> lifetimeGroupPars() const;

    // Physical fractions computed from the parameter vector.
    struct Fractions {
        std::vector<double> resolution;   // w_i, sum = 1
        std::vector<double> lifetime;     // I_j, p-Ps entry = total p-Ps intensity
    };
    Fractions fractions(const double* p) const;

    // Initial values of the fraction/intensity parameters from the components.
    void setInitialIntensityParameters(std::vector<double>& p) const;

    std::vector<std::string> parameterNames() const;

    // Restrict the evaluation to a single lifetime component (for drawing).
    // -1 = all components.
    void setComponentMask(int j) { componentMask_ = j; }
    void setIncludeBackground(bool b) { includeBackground_ = b; }

    double operator()(const double* x, const double* p) const;

private:
    void computeFractions(const double* p, Fractions& f) const;

    std::vector<ResolutionComponent> resolution_;
    std::vector<LifetimeComponent> lifetimes_;
    double binWidth_;
    double oPsThreshold_;

    ResolutionLayout resLayout_;
    std::vector<int> lifeFree_, lifeAbs_;
    int pPs_ = -1;

    std::vector<int> widthPar_;
    int nWidths_ = 0;

    int componentMask_ = -1;
    bool includeBackground_ = true;

    // Fractions are computed once per parameter set, not once per bin.
    // (Not thread-safe: do not use with multithreaded TF1 evaluation.)
    mutable std::vector<double> cachedPars_;
    mutable Fractions cached_;
};

} // namespace pals

#endif
