#ifndef PALS_SPECTRUM_H
#define PALS_SPECTRUM_H

#include "Config.h"

#include <cstddef>
#include <iosfwd>
#include <string>
#include <vector>

namespace pals {

// Binned lifetime spectrum with uniform bins (ROOT-independent).
struct Spectrum {
    double min = 0.0;            // low edge of bin 0 [ns]
    double binWidth = 0.0;       // [ns]
    std::vector<double> counts;

    std::size_t size() const { return counts.size(); }
    double center(std::size_t i) const { return min + (i + 0.5) * binWidth; }
    double max() const { return min + size() * binWidth; }
};

// Reads one whitespace-separated column of a text file. Lines starting with
// '#' and non-numeric / non-finite values (nan, inf) are skipped.
std::vector<double> readTextColumn(const std::string& path, int column);

// Histogram of time differences with the binning from the config.
Spectrum binTimes(const std::vector<double>& times, const HistogramConfig& h);

// Text histogram: one value (counts) per bin, first bin starting at h.min.
Spectrum spectrumFromCounts(std::vector<double> counts, const HistogramConfig& h);

// Preliminary estimates used to set the fit range and initial parameters.
// Bin indices are 0-based indices of Spectrum::counts.
struct Estimates {
    std::size_t peakBin = 0;
    std::size_t firstBin = 0;        // first bin above min_counts_first_bin
    std::size_t lastBin = 0;         // last bin above min_counts_last_bin
    std::size_t backgroundFrom = 0, backgroundTo = 0;   // inclusive
    double background = 0.0;         // mean counts per bin
    double backgroundError = 0.0;    // standard error of the mean
    std::size_t fitFrom = 0, fitTo = 0;                 // inclusive
    double peakPosition = 0.0;       // [ns]
    double peakCounts = 0.0;
    double area = 0.0;               // counts above background, firstBin..lastBin
    double sigma = 0.0;              // from the left half-width of the peak [ns]
};

Estimates estimate(const Spectrum& s, const EstimationConfig& cfg);
void printEstimates(std::ostream& os, const Spectrum& s, const Estimates& e);

} // namespace pals

#endif
