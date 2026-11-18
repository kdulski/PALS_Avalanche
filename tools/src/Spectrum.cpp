#include "Spectrum.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace pals {

std::vector<double> readTextColumn(const std::string& path, int column)
{
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Cannot open data file " + path);

    std::vector<double> values;
    std::size_t skipped = 0;
    std::string line;
    while (std::getline(in, line)) {
        const char* c = line.c_str();
        while (*c == ' ' || *c == '\t') ++c;
        if (*c == '\0' || *c == '#') continue;
        // jump to the requested column
        for (int col = 0; col < column && *c; ++col) {
            while (*c && *c != ' ' && *c != '\t') ++c;
            while (*c == ' ' || *c == '\t') ++c;
        }
        char* end = nullptr;
        const double v = std::strtod(c, &end);
        if (end == c || !std::isfinite(v)) { ++skipped; continue; }
        values.push_back(v);
    }
    std::cout << "Read " << values.size() << " values from column " << column << " of " << path;
    if (skipped) std::cout << " (" << skipped << " lines skipped: missing column, nan or inf)";
    std::cout << "\n";
    if (values.empty()) throw std::runtime_error("No values read from " + path);
    return values;
}

Spectrum binTimes(const std::vector<double>& times, const HistogramConfig& h)
{
    Spectrum s;
    s.min = h.min;
    s.binWidth = h.binWidth;
    const double nExact = (h.max - h.min) / h.binWidth;
    const auto n = static_cast<std::size_t>(std::llround(nExact));
    if (n == 0) throw std::runtime_error("Histogram range is smaller than one bin");
    if (std::fabs(nExact - n) > 1e-6 * nExact)
        std::cerr << "Warning: histogram range is not a multiple of bin_width; using " << n
                  << " bins, upper edge " << s.min + n * s.binWidth << " ns\n";
    s.counts.assign(n, 0.0);
    std::size_t outside = 0;
    for (double t : times) {
        const double x = (t - s.min) / s.binWidth;
        if (x < 0.0 || x >= static_cast<double>(n)) { ++outside; continue; }
        s.counts[static_cast<std::size_t>(x)] += 1.0;
    }
    std::cout << "Histogram: " << n << " bins of " << s.binWidth << " ns in [" << s.min << ", " << s.max()
              << "] ns, " << times.size() - outside << " entries, " << outside << " outside the range\n";
    return s;
}

Spectrum spectrumFromCounts(std::vector<double> counts, const HistogramConfig& h)
{
    Spectrum s;
    s.min = h.min;
    s.binWidth = h.binWidth;
    s.counts = std::move(counts);
    std::cout << "Histogram: " << s.size() << " bins of " << s.binWidth << " ns starting at " << s.min << " ns\n";
    return s;
}

Estimates estimate(const Spectrum& s, const EstimationConfig& cfg)
{
    const auto& c = s.counts;
    const std::size_t n = s.size();
    if (n < 3) throw std::runtime_error("Spectrum has too few bins");
    Estimates e;

    e.peakBin = static_cast<std::size_t>(std::max_element(c.begin(), c.end()) - c.begin());
    e.peakPosition = s.center(e.peakBin);
    e.peakCounts = c[e.peakBin];

    auto first = std::find_if(c.begin(), c.end(), [&](double v) { return v > cfg.minCountsFirstBin; });
    auto last = std::find_if(c.rbegin(), c.rend(), [&](double v) { return v > cfg.minCountsLastBin; });
    if (first == c.end() || last == c.rend())
        throw std::runtime_error("No bin above min_counts_first_bin / min_counts_last_bin");
    e.firstBin = static_cast<std::size_t>(first - c.begin());
    e.lastBin = n - 1 - static_cast<std::size_t>(last - c.rbegin());

    // --- background: mean of nBins bins, skipBins away from the end of the spectrum,
    //     or the bins with centres in an explicit window
    const auto& bg = cfg.background;
    long from, to;
    if (bg.useRange) {
        from = static_cast<long>(std::ceil((bg.from - s.min) / s.binWidth - 0.5));
        to = static_cast<long>(std::floor((bg.to - s.min) / s.binWidth - 0.5));
        from = std::max(from, 0L);
        to = std::min(to, static_cast<long>(n) - 1);
    } else {
        const long nb = bg.nBins, skip = bg.skipBins;
        if (bg.rightSide) { to = static_cast<long>(e.lastBin) - skip - 1; from = to - nb + 1; }
        else              { from = static_cast<long>(e.firstBin) + skip + 1; to = from + nb - 1; }
    }
    if (from < 0 || to >= static_cast<long>(n) || to - from < 1)
        throw std::runtime_error("Background window does not fit into the spectrum (or has < 2 bins)");
    const long nb = to - from + 1;
    e.backgroundFrom = static_cast<std::size_t>(from);
    e.backgroundTo = static_cast<std::size_t>(to);
    double sum = 0.0;
    for (long i = from; i <= to; ++i) sum += c[i];
    e.background = sum / nb;
    double ss = 0.0;
    for (long i = from; i <= to; ++i) ss += (c[i] - e.background) * (c[i] - e.background);
    e.backgroundError = std::sqrt(ss / (nb * (nb - 1.0)));

    const double height = e.peakCounts - e.background;
    if (height <= 0.0) throw std::runtime_error("Peak is not above the background");

    // --- fit range: start where the leading edge drops below fitStartFraction of the peak,
    //     end at the first bin with centre >= fitEnd (or at the last bin of the spectrum)
    e.fitFrom = 0;
    if (cfg.useFitStart) {
        while (e.fitFrom + 1 < n && s.center(e.fitFrom) < cfg.fitStart) ++e.fitFrom;
    } else {
        for (std::size_t i = e.peakBin; i-- > 0;) {
            if (c[i] - e.background < cfg.fitStartFraction * height) { e.fitFrom = i; break; }
        }
    }
    e.fitTo = e.lastBin;
    for (std::size_t i = e.peakBin; i <= e.lastBin; ++i) {
        if (s.center(i) >= cfg.fitEnd) { e.fitTo = i; break; }
    }
    if (e.fitTo <= e.peakBin) throw std::runtime_error("fit_end lies before the peak of the spectrum");
    if (e.fitFrom >= e.fitTo) throw std::runtime_error("fit start is not before fit_end");

    // --- area: background-subtracted counts of the whole spectrum
    double area = 0.0;
    for (std::size_t i = e.firstBin; i <= e.lastBin; ++i) area += c[i] - e.background;
    e.area = area;

    // --- sigma of the prompt peak from its left half-width (Gaussian: HWHM = sigma*sqrt(2 ln 2))
    e.sigma = s.binWidth;
    for (std::size_t i = e.peakBin; i-- > 0;) {
        const double yi = c[i] - e.background;
        if (yi < 0.5 * height) {
            const double yn = c[i + 1] - e.background;
            const double xHalf = s.center(i) + (yn > yi ? (0.5 * height - yi) / (yn - yi) : 0.0) * s.binWidth;
            e.sigma = std::max(s.binWidth, (e.peakPosition - xHalf) / std::sqrt(2.0 * std::log(2.0)));
            break;
        }
    }
    return e;
}

void printEstimates(std::ostream& os, const Spectrum& s, const Estimates& e)
{
    auto bin = [&](std::size_t i) {
        std::ostringstream o;
        o << "bin " << i << " (" << std::setprecision(5) << s.center(i) << " ns)";
        return o.str();
    };
    os << "Peak:              " << bin(e.peakBin) << ", " << e.peakCounts << " counts\n"
       << "Spectrum:          " << bin(e.firstBin) << " .. " << bin(e.lastBin) << "\n"
       << "Background:        " << e.background << " +- " << e.backgroundError << " counts/bin, from "
       << bin(e.backgroundFrom) << " .. " << bin(e.backgroundTo) << "\n"
       << "Fit range:         " << bin(e.fitFrom) << " .. " << bin(e.fitTo) << "\n"
       << "Area:              " << e.area << " counts above background\n"
       << "Sigma (left HWHM): " << e.sigma << " ns\n";
}

} // namespace pals
