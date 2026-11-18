#ifndef PALS_DISTRIBUTION_H
#define PALS_DISTRIBUTION_H

// Lifetime distribution on a grid, after the discrete fit (ROOT-independent):
//   - NNLS (Lawson-Hanson) averaged over parametric bootstrap replicas ("bagged")
//   - MaxEnt with a flat default model
// and a comparison of the two (features determined by the data vs by the method).
//
// Both use the resolution, t0 and background of the discrete fit, the same fit range
// and Pearson weights 1/model. MaxEnt's alpha is chosen so that
// chi2 = chi2(NNLS) + sqrt(2 * bins) (discrepancy principle).

#include "Config.h"

#include <iosfwd>
#include <string>
#include <vector>

namespace pals {

struct DistributionProblem {
    std::vector<double> t, y, weight;    // fit-range bins: centre, counts, 1/max(model, 1)
    double binWidth = 0.0;
    struct Gaussian { double weight, sigma, t0; };
    std::vector<Gaussian> resolution;
    double background = 0.0;             // from the discrete fit (start value)
};

// Summary of a distribution in a lifetime region [low, high]
struct RegionSummary {
    double low = 0.0, high = 0.0;
    double intensity = 0.0;              // fraction of all components (background excluded)
    double centroid = 0.0;               // intensity-weighted mean tau [ns]
    double width = 0.0;                  // intensity-weighted sd of tau [ns]
    int peaks = 0;                       // local maxima of the density per ln(tau) above 20 % of the maximum
};

// Spread of a summary over bootstrap replicas (bagged NNLS)
struct RegionSpread {
    double intensitySd = 0.0, centroidSd = 0.0;
    double widthQ16 = 0.0, widthQ84 = 0.0;
    double multiPeakFraction = 0.0;      // replicas with >= 2 peaks
};

struct DistributionEstimate {
    std::string method;
    std::vector<double> x;               // counts per grid point
    std::vector<double> densityLow, densityHigh;   // 16 / 84 % bootstrap band of the density (bagged NNLS)
    double background = 0.0;
    double chi2 = 0.0;
    double alpha = 0.0;                  // MaxEnt regularisation
    std::vector<RegionSummary> regions;
    std::vector<RegionSpread> spreads;   // bagged NNLS only
};

struct DistributionAnalysis {
    std::vector<double> grid;            // lifetimes [ns], log-spaced
    double chi2Nnls = 0.0;               // unregularised NNLS (reference for alpha)
    int bins = 0;
    std::vector<DistributionEstimate> estimates;
    std::vector<std::string> comparison; // text: consensus and flags
};

std::vector<double> logGrid(double lo, double hi, int n);

// density per ln(tau), normalised to unit total intensity
std::vector<double> densityPerLnTau(const std::vector<double>& grid, const std::vector<double>& x);

RegionSummary summarize(const std::vector<double>& grid, const std::vector<double>& x, double low, double high);

DistributionAnalysis analyzeDistribution(const DistributionProblem& problem, const DistributionConfig& cfg,
                                         const std::vector<std::pair<double, double>>& regions);

void writeDistributionReport(std::ostream& os, const DistributionAnalysis& a, const DistributionConfig& cfg,
                             double counts);

} // namespace pals

#endif
