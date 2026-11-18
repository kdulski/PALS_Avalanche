#include "Distribution.h"
#include "Emg.h"
#include "Nnls.h"
#include "TextReport.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <ostream>
#include <random>
#include <sstream>

namespace pals {

std::vector<double> logGrid(double lo, double hi, int n)
{
    std::vector<double> g(n);
    for (int k = 0; k < n; ++k) g[k] = std::exp(std::log(lo) + (std::log(hi) - std::log(lo)) * k / (n - 1));
    return g;
}

std::vector<double> densityPerLnTau(const std::vector<double>& grid, const std::vector<double>& x)
{
    const double total = std::accumulate(x.begin(), x.end(), 0.0);
    const double dl = std::log(grid[1] / grid[0]);
    std::vector<double> d(x.size(), 0.0);
    if (total > 0.0)
        for (std::size_t k = 0; k < x.size(); ++k) d[k] = x[k] / total / dl;
    return d;
}

RegionSummary summarize(const std::vector<double>& grid, const std::vector<double>& x, double low, double high)
{
    RegionSummary s;
    s.low = low;
    s.high = high;
    const double total = std::accumulate(x.begin(), x.end(), 0.0);
    double sum = 0.0, sumT = 0.0, dmax = 0.0;
    std::vector<std::size_t> idx;
    for (std::size_t k = 0; k < grid.size(); ++k) {
        if (grid[k] < low || grid[k] > high) continue;
        idx.push_back(k);
        sum += x[k];
        sumT += x[k] * grid[k];
        dmax = std::max(dmax, x[k]);
    }
    if (total <= 0.0 || sum <= 0.0) return s;
    s.intensity = sum / total;
    s.centroid = sumT / sum;
    double var = 0.0;
    for (std::size_t k : idx) var += x[k] * (grid[k] - s.centroid) * (grid[k] - s.centroid);
    s.width = std::sqrt(var / sum);
    // uniform grid in ln(tau): density ~ x
    for (std::size_t q = 1; q + 1 < idx.size(); ++q) {
        const double v = x[idx[q]];
        if (v >= x[idx[q - 1]] && v > x[idx[q + 1]] && v > 0.2 * dmax) ++s.peaks;
    }
    return s;
}

namespace {

struct Kernel {
    ColMatrix Kb;                        // [K, 1] unweighted, bins x (n + 1)
    ColMatrix A;                         // sqrt(w) * [K, 1]
    std::vector<double> sw;              // sqrt(weights)
};

Kernel buildKernel(const DistributionProblem& pr, const std::vector<double>& grid)
{
    const std::size_t m = pr.t.size(), n = grid.size();
    Kernel k;
    k.Kb = ColMatrix(m, n + 1);
    k.A = ColMatrix(m, n + 1);
    k.sw.resize(m);
    for (std::size_t i = 0; i < m; ++i) k.sw[i] = std::sqrt(pr.weight[i]);
    for (std::size_t j = 0; j < n; ++j)
        for (std::size_t i = 0; i < m; ++i) {
            double v = 0.0;
            for (const auto& g : pr.resolution) v += g.weight * pals_emg(pr.t[i], g.sigma, g.t0, grid[j]);
            k.Kb(i, j) = pr.binWidth * v;
        }
    for (std::size_t i = 0; i < m; ++i) k.Kb(i, n) = 1.0;
    for (std::size_t j = 0; j <= n; ++j)
        for (std::size_t i = 0; i < m; ++i) k.A(i, j) = k.sw[i] * k.Kb(i, j);
    return k;
}

double weightedChi2(const Kernel& k, const DistributionProblem& pr, const std::vector<double>& v)
{
    double chi2 = 0.0;
    for (std::size_t i = 0; i < pr.t.size(); ++i) {
        double m = 0.0;
        for (std::size_t j = 0; j < v.size(); ++j) m += k.Kb(i, j) * v[j];
        chi2 += pr.weight[i] * (pr.y[i] - m) * (pr.y[i] - m);
    }
    return chi2;
}

double quantile(std::vector<double> v, double q)
{
    std::sort(v.begin(), v.end());
    const double pos = q * (v.size() - 1);
    const std::size_t i = static_cast<std::size_t>(pos);
    return i + 1 < v.size() ? v[i] + (pos - i) * (v[i + 1] - v[i]) : v.back();
}

double sd(const std::vector<double>& v)
{
    const double m = std::accumulate(v.begin(), v.end(), 0.0) / v.size();
    double s = 0.0;
    for (double x : v) s += (x - m) * (x - m);
    return std::sqrt(s / std::max<std::size_t>(1, v.size() - 1));
}

// MaxEnt: minimise chi2/2 - alpha * S(x), S = sum(x - m - x ln(x/m)), x > 0, background free;
// damped Newton (Hessian Q + alpha diag(1/x) is positive definite).
struct MaxEntSolver {
    const Kernel& k;
    const DistributionProblem& pr;
    std::size_t n;
    std::vector<double> Q, c, m;
    double d = 0.0;

    MaxEntSolver(const Kernel& kk, const DistributionProblem& p, std::vector<double> defaultModel)
        : k(kk), pr(p), n(defaultModel.size()), m(std::move(defaultModel))
    {
        const std::size_t N = n + 1, rows = pr.t.size();
        Q.assign(N * N, 0.0);
        c.assign(N, 0.0);
        for (std::size_t a = 0; a < N; ++a) {
            const double* ca = k.Kb.col(a);
            for (std::size_t b = a; b < N; ++b) {
                const double* cb = k.Kb.col(b);
                double s = 0.0;
                for (std::size_t i = 0; i < rows; ++i) s += pr.weight[i] * ca[i] * cb[i];
                Q[a * N + b] = Q[b * N + a] = s;
            }
            double s = 0.0;
            for (std::size_t i = 0; i < rows; ++i) s += pr.weight[i] * ca[i] * pr.y[i];
            c[a] = s;
        }
        for (std::size_t i = 0; i < rows; ++i) d += pr.weight[i] * pr.y[i] * pr.y[i];
    }

    double chi2(const std::vector<double>& v) const
    {
        const std::size_t N = n + 1;
        double s = d;
        for (std::size_t a = 0; a < N; ++a) {
            double qa = 0.0;
            for (std::size_t b = 0; b < N; ++b) qa += Q[a * N + b] * v[b];
            s += v[a] * qa - 2.0 * c[a] * v[a];
        }
        return s;
    }

    double objective(const std::vector<double>& v, double alpha) const
    {
        double S = 0.0;
        for (std::size_t j = 0; j < n; ++j) S += v[j] - m[j] - v[j] * std::log(v[j] / m[j]);
        return 0.5 * chi2(v) - alpha * S;
    }

    std::vector<double> solve(double alpha, std::vector<double> v) const
    {
        const std::size_t N = n + 1;
        for (int it = 0; it < 300; ++it) {
            std::vector<double> g(N), H = Q, step;
            for (std::size_t a = 0; a < N; ++a) {
                double qa = 0.0;
                for (std::size_t b = 0; b < N; ++b) qa += Q[a * N + b] * v[b];
                g[a] = qa - c[a];
            }
            for (std::size_t j = 0; j < n; ++j) {
                g[j] += alpha * std::log(v[j] / m[j]);
                H[j * N + j] += alpha / v[j];
            }
            std::vector<double> minus(N);
            for (std::size_t a = 0; a < N; ++a) minus[a] = -g[a];
            if (!choleskySolve(H, step, minus, N)) break;
            double gd = 0.0;
            for (std::size_t a = 0; a < N; ++a) gd += g[a] * step[a];
            if (-gd < 1e-10 * std::max(1.0, std::fabs(objective(v, alpha)))) break;
            double s = 1.0;
            for (std::size_t j = 0; j < n; ++j)
                if (step[j] < 0.0) s = std::min(s, 0.9 * v[j] / -step[j]);
            const double f0 = objective(v, alpha);
            std::vector<double> trial(N);
            for (int ls = 0; ls < 50; ++ls) {
                for (std::size_t a = 0; a < N; ++a) trial[a] = v[a] + s * step[a];
                if (objective(trial, alpha) <= f0 + 1e-4 * s * gd) break;
                s *= 0.5;
            }
            v = trial;
        }
        return v;
    }
};

} // namespace

DistributionAnalysis analyzeDistribution(const DistributionProblem& pr, const DistributionConfig& cfg,
                                         const std::vector<std::pair<double, double>>& regions)
{
    DistributionAnalysis out;
    out.grid = logGrid(cfg.gridMin, cfg.gridMax, cfg.gridPoints);
    out.bins = static_cast<int>(pr.t.size());
    const std::size_t n = out.grid.size(), rows = pr.t.size();
    const Kernel k = buildKernel(pr, out.grid);

    std::vector<double> b(rows);
    for (std::size_t i = 0; i < rows; ++i) b[i] = k.sw[i] * pr.y[i];
    const NnlsResult single = nnls(k.A, b);
    out.chi2Nnls = single.residualNorm * single.residualNorm;
    const double target = out.chi2Nnls + std::sqrt(2.0 * rows);
    auto summaries = [&](const std::vector<double>& x) {
        std::vector<RegionSummary> s;
        for (const auto& r : regions) s.push_back(summarize(out.grid, x, r.first, r.second));
        return s;
    };

    if (cfg.nnlsBagged) {
        std::cout << "Distribution: NNLS (Lawson-Hanson) with " << cfg.bootstrap << " bootstrap replicas\n";
        // parametric bootstrap from the spectrum fitted by NNLS itself
        std::vector<double> model(rows, 0.0);
        for (std::size_t j = 0; j <= n; ++j)
            for (std::size_t i = 0; i < rows; ++i) model[i] += k.Kb(i, j) * single.x[j];
        std::mt19937_64 rng(cfg.seed);
        DistributionEstimate e;
        e.method = "NNLS bagged";
        e.x.assign(n, 0.0);
        std::vector<std::vector<double>> densities;
        std::vector<std::vector<RegionSummary>> reps;
        for (int r = 0; r < cfg.bootstrap; ++r) {
            std::vector<double> yb(rows);
            for (std::size_t i = 0; i < rows; ++i) {
                std::poisson_distribution<long long> pois(std::max(model[i], 0.0));
                yb[i] = k.sw[i] * static_cast<double>(pois(rng));
            }
            const NnlsResult rb = nnls(k.A, yb);
            std::vector<double> xb(rb.x.begin(), rb.x.begin() + n);
            for (std::size_t j = 0; j < n; ++j) e.x[j] += xb[j] / cfg.bootstrap;
            e.background += rb.x[n] / cfg.bootstrap;
            densities.push_back(densityPerLnTau(out.grid, xb));
            reps.push_back(summaries(xb));
        }
        e.densityLow.resize(n);
        e.densityHigh.resize(n);
        for (std::size_t j = 0; j < n; ++j) {
            std::vector<double> v;
            for (const auto& d : densities) v.push_back(d[j]);
            e.densityLow[j] = quantile(v, 0.16);
            e.densityHigh[j] = quantile(v, 0.84);
        }
        std::vector<double> v(e.x);
        v.push_back(e.background);
        e.chi2 = weightedChi2(k, pr, v);
        e.regions = summaries(e.x);
        for (std::size_t r = 0; r < regions.size(); ++r) {
            std::vector<double> I, C, W;
            int multi = 0;
            for (const auto& rep : reps) {
                I.push_back(rep[r].intensity);
                C.push_back(rep[r].centroid);
                W.push_back(rep[r].width);
                multi += rep[r].peaks >= 2;
            }
            RegionSpread s;
            s.intensitySd = sd(I);
            s.centroidSd = sd(C);
            s.widthQ16 = quantile(W, 0.16);
            s.widthQ84 = quantile(W, 0.84);
            s.multiPeakFraction = static_cast<double>(multi) / reps.size();
            e.spreads.push_back(s);
        }
        out.estimates.push_back(e);
    }

    if (cfg.maxent) {
        std::cout << "Distribution: MaxEnt (flat default model, alpha from chi2 = chi2_NNLS + sqrt(2 bins))\n";
        double total = 0.0;
        for (std::size_t j = 0; j < n; ++j) total += single.x[j];
        const MaxEntSolver solver(k, pr, std::vector<double>(n, std::max(total, 1.0) / n));
        std::vector<double> start(solver.m);
        start.push_back(single.x[n]);
        double lo = -8.0, hi = 3.0;
        std::vector<double> best;
        double bestAlpha = 0.0;
        std::vector<double> warm = start;
        for (int it = 0; it < 16; ++it) {
            const double mid = 0.5 * (lo + hi);
            const std::vector<double> v = solver.solve(std::pow(10.0, mid), warm);
            if (solver.chi2(v) <= target) { lo = mid; best = v; bestAlpha = std::pow(10.0, mid); warm = v; }
            else hi = mid;
        }
        if (best.empty()) { bestAlpha = std::pow(10.0, lo); best = solver.solve(bestAlpha, warm); }
        DistributionEstimate e;
        e.method = "MaxEnt";
        e.x.assign(best.begin(), best.begin() + n);
        e.background = best[n];
        e.chi2 = solver.chi2(best);
        e.alpha = bestAlpha;
        e.regions = summaries(e.x);
        out.estimates.push_back(e);
    }

    // --- comparison: what the two methods (sparse vs smooth prior) agree on
    const DistributionEstimate* bag = nullptr;
    const DistributionEstimate* me = nullptr;
    for (const auto& e : out.estimates) (e.method == "MaxEnt" ? me : bag) = &e;
    if (bag && me) {
        const double step = std::log(out.grid[1] / out.grid[0]);
        for (std::size_t r = 0; r < regions.size(); ++r) {
            const auto& B = bag->regions[r];
            const auto& M = me->regions[r];
            const auto& S = bag->spreads[r];
            std::ostringstream o;
            o << std::setprecision(4);
            o << "Region " << B.low << " - " << B.high << " ns:";
            out.comparison.push_back(o.str());
            std::vector<std::string> robust, open;
            auto line = [&](const std::string& s) { out.comparison.push_back("  " + s); };
            {
                std::ostringstream l;
                const bool ok = std::fabs(B.intensity - M.intensity) <= std::max(2.0 * S.intensitySd, 0.002);
                l << std::setprecision(4) << "intensity: NNLS-bagged " << 100 * B.intensity << " +- " << 100 * S.intensitySd
                  << " %, MaxEnt " << 100 * M.intensity << " %  -> " << (ok ? "consistent" : "DIFFERENT");
                line(l.str());
                (ok ? robust : open).push_back("intensity");
            }
            {
                std::ostringstream l;
                const double tol = 2.0 * S.centroidSd + B.centroid * step;
                const bool ok = std::fabs(B.centroid - M.centroid) <= tol;
                l << std::setprecision(4) << "centroid:  NNLS-bagged " << B.centroid << " +- " << S.centroidSd << " ns, MaxEnt "
                  << M.centroid << " ns  -> " << (ok ? "consistent" : "DIFFERENT");
                line(l.str());
                (ok ? robust : open).push_back("centroid");
            }
            {
                std::ostringstream l;
                l << std::setprecision(3) << "width (sd of tau): NNLS-bagged " << B.width << " [" << S.widthQ16 << ", "
                  << S.widthQ84 << "] ns, MaxEnt " << M.width << " ns  -> ";
                if (M.width > S.widthQ84) {
                    l << "MaxEnt broader: the width is set mainly by the regularisation; the data are also "
                         "compatible with a narrower distribution";
                    open.push_back("width");
                } else if (M.width < S.widthQ16) {
                    l << "MaxEnt narrower than the bootstrap NNLS range";
                    open.push_back("width");
                } else {
                    l << "consistent";
                    robust.push_back("width");
                }
                line(l.str());
            }
            {
                std::ostringstream l;
                l << std::setprecision(2) << "maxima: NNLS-bagged >= 2 maxima in " << 100 * S.multiPeakFraction
                  << " % of the bootstrap replicas, MaxEnt " << M.peaks << "  -> ";
                const bool nnlsMulti = S.multiPeakFraction >= 0.8, nnlsSingle = S.multiPeakFraction <= 0.2;
                const bool ok = (nnlsMulti && M.peaks >= 2) || (nnlsSingle && M.peaks <= 1);
                l << (ok ? "consistent" : "number of maxima NOT determined by the data");
                line(l.str());
                (ok ? robust : open).push_back("number of maxima");
            }
            auto join = [](const std::vector<std::string>& v) {
                std::string s;
                for (std::size_t q = 0; q < v.size(); ++q) s += (q ? ", " : "") + v[q];
                return s.empty() ? std::string("-") : s;
            };
            line("same with sparse and smooth prior:  " + join(robust));
            line("depends on the prior (not determined by the data): " + join(open));
        }
    }
    return out;
}

void writeDistributionReport(std::ostream& os, const DistributionAnalysis& a, const DistributionConfig& cfg,
                             double counts)
{
    os << "\n=== Lifetime distribution (grid " << cfg.gridMin << " - " << cfg.gridMax << " ns, " << cfg.gridPoints
       << " points, log-spaced) ===\n"
       << statisticsNote(counts, cfg.countsLow, cfg.countsGood, "lifetime distribution")
       << "Fit range " << a.bins << " bins; chi2 of unregularised NNLS (Pearson weights): " << a.chi2Nnls << "\n";
    for (const auto& e : a.estimates) {
        os << "\n" << e.method << ": chi2 " << e.chi2 << ", background " << e.background;
        if (e.method == "MaxEnt") os << ", alpha " << e.alpha;
        if (e.method == "NNLS bagged") os << ", " << cfg.bootstrap << " bootstrap replicas";
        os << "\n  " << std::left << std::setw(16) << "region [ns]" << std::setw(18) << "intensity [%]"
           << std::setw(18) << "centroid [ns]" << std::setw(26) << "width (sd) [ns]" << "maxima\n";
        for (std::size_t r = 0; r < e.regions.size(); ++r) {
            const auto& s = e.regions[r];
            std::ostringstream reg, inten, cen, wid;
            reg << std::setprecision(3) << s.low << " - " << s.high;
            inten << std::setprecision(4) << 100 * s.intensity;
            cen << std::setprecision(4) << s.centroid;
            wid << std::setprecision(3) << s.width;
            if (!e.spreads.empty()) {
                inten << " +- " << std::setprecision(2) << 100 * e.spreads[r].intensitySd;
                cen << " +- " << std::setprecision(2) << e.spreads[r].centroidSd;
                wid << " [" << std::setprecision(3) << e.spreads[r].widthQ16 << ", " << e.spreads[r].widthQ84 << "]";
            }
            os << "  " << std::setw(16) << reg.str() << std::setw(18) << inten.str() << std::setw(18) << cen.str()
               << std::setw(26) << wid.str() << s.peaks << "\n";
        }
        os << std::right;
    }
    if (!a.comparison.empty()) {
        os << "\nComparison NNLS-bagged (sparse prior) vs MaxEnt (smooth prior):\n";
        for (const auto& l : a.comparison) os << l << "\n";
        os << "Features on which the two disagree are not determined by the spectrum. Agreement is necessary but\n"
              "not sufficient: e.g. two close components and a continuous distribution with the same mean and width\n"
              "give the same result. Report shapes only together with a simulation of the measurement (same\n"
              "resolution and statistics).\n";
    }
    os << "The distributions (density per ln tau) are stored in the ROOT file.\n";
}

} // namespace pals
