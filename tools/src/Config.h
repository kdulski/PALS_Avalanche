#ifndef PALS_CONFIG_H
#define PALS_CONFIG_H

#include "Components.h"

#include <string>
#include <vector>

namespace pals {

enum class InputFormat { Times, Histogram, Root };

struct InputConfig {
    InputFormat format = InputFormat::Times;
    int column = 0;              // text files: 0-based column with the value
    std::string rootObject;      // ROOT files: path of the histogram, e.g. "dir/histo"
};

struct HistogramConfig {
    double binWidth = 0.0;       // [ns]
    double min = 0.0;            // low edge of the first bin [ns]
    double max = 0.0;            // high edge of the last bin [ns] (times only)
};

struct BackgroundConfig {
    bool rightSide = true;       // estimate from the right (long-time) end of the spectrum
    int nBins = 100;             // number of bins averaged
    int skipBins = 30;           // bins skipped at the very end of the spectrum
    double nSigma = 3.0;         // background may vary by +- nSigma * (its uncertainty)
    bool useRange = false;       // explicit window [from, to] in ns instead of side/bins/skip
    double from = 0.0, to = 0.0;
};

struct EstimationConfig {
    double minCountsFirstBin = 8;  // spectrum starts at the first bin above this
    double minCountsLastBin = 8;   // spectrum ends at the last bin above this
    BackgroundConfig background;
    double fitStartFraction = 0.3; // fit starts where (N - B) drops below this fraction of the peak
    bool useFitStart = false;      // ... or at an explicit time [ns]
    double fitStart = 0.0;
    double fitEnd = 40.0;          // [ns]
};

struct FitConfig {
    int iterations = 3;
    double partialVariation = 0.05;     // fraction (config gives percent)
    double tauMin = 0.08, tauMax = 142.0;   // limits of free lifetimes [ns]
    double sigmaMin = 0.01, sigmaMax = 5.0; // limits of free sigmas [ns]
    double offsetWindow = 5.0;          // t0 in [peak - k*sigma, peak + k*sigma]
    double oPsThreshold = 0.7;          // components with tau above are o-Ps (for p-Ps) [ns]
    std::string firstPassOptions = "RMW";
    std::string options = "RM";
    double widthMin = 0.0, widthMax = 1.0;   // limits of component widths (sd of ln tau)
    bool profileWidths = true;          // profile-likelihood intervals for free widths
};

// Optional analysis of the lifetime distribution after the discrete fit:
// NNLS (Lawson-Hanson) averaged over bootstrap replicas and MaxEnt on a grid
// of lifetimes, plus a comparison of the two.
struct DistributionConfig {
    bool enabled = false;
    bool nnlsBagged = true;
    bool maxent = true;
    double gridMin = 0.08, gridMax = 15.0;   // [ns], log-spaced grid
    int gridPoints = 160;
    int bootstrap = 30;
    unsigned seed = 12345;
    std::vector<std::pair<double, double>> regions;   // summaries; default: [o-Ps threshold, grid max]
    double countsLow = 5e6, countsGood = 1e7;         // statistics thresholds for the warning
};

enum class ModelType { Lifetimes, Formula };

// Parameter of a formula model. The initial value is either a number or a
// preliminary estimate ("estimate:<name>[*factor]", see ModelSetup.h).
struct FormulaParameter {
    std::string name;
    double value = 0.0;
    std::string estimate;        // empty = use value
    double estimateFactor = 1.0;
    bool hasLimits = false;
    double lower = 0.0, upper = 0.0;
    ParMode mode = ParMode::Free;
    int stage = 0;               // first fit stage in which the parameter is released
};

struct ModelConfig {
    ModelType type = ModelType::Lifetimes;
    std::string formula;                  // TFormula expression in x and [named] parameters
    bool convolve = true;                 // convolve the formula with the resolution function
    std::vector<FormulaParameter> parameters;
    double pointsPerSigma = 5.0;          // grid step = min(sigma)/pointsPerSigma ...
    double step = 0.0;                    // ... unless given explicitly [ns]
    double windowSigma = 8.0;             // Gaussian truncated at +- windowSigma * sigma
};

struct OutputConfig {
    std::string directory = "Results";
    bool timestamp = true;
};

struct Config {
    InputConfig input;
    HistogramConfig histogram;
    EstimationConfig estimation;
    std::vector<ResolutionComponent> resolution;
    std::vector<LifetimeComponent> lifetimes;
    ModelConfig model;
    FitConfig fit;
    DistributionConfig distribution;
    OutputConfig output;

    std::string description;
    std::string sourcePath;
    std::string sourceText;      // file content, stored in the output for provenance

    static Config load(const std::string& path);
    static Config parse(const std::string& jsonText, const std::string& sourcePath = "");

    // JSON text of the configuration (all settings written explicitly); parse(toJson())
    // gives the same configuration.
    std::string toJson() const;
    void save(const std::string& path) const;
};

} // namespace pals

#endif
