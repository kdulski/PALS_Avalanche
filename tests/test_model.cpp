// Tests of the ROOT-independent part:
//  1. hyperspherical parametrisation (round trip, comparison with the original formulas)
//  2. full model vs the original DiscreteFitFunctionNoPs / DiscreteFitFunctionPs
//  3. EMG accuracy and normalisation
//  4. benchmark: original model vs the new one
//  5. formula model (numerical convolution / no convolution) vs the lifetime model
//  6. component width (log-normal), Gauss-Hermite quadrature, NNLS, distribution analysis
//     on simulated spectra (S0: one o-Ps component, S1: 2.0 + 2.5 ns, S2: log-normal)
//  7. (optional) preliminary estimates, formula configs, config round trip (toJson -> parse)
//     and explicit fit start / background window on a real spectrum:
//       test_model <data file> <config.json>

#include "Config.h"
#include "Distribution.h"
#include "Nnls.h"
#include "Quadrature.h"
#include "Convolution.h"
#include "Emg.h"
#include "FormulaModel.h"
#include "ModelSetup.h"
#include "Intensities.h"
#include "LifetimeModel.h"
#include "Spectrum.h"
#include "reference_original.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <random>
#include <vector>

using namespace pals;

static int failures = 0;
static void check(bool ok, const char* what, double value)
{
    std::printf("  %-62s %-4s (%.3g)\n", what, ok ? "OK" : "FAIL", value);
    if (!ok) ++failures;
}

// ---------------------------------------------------------------------------
static void testHyperspherical()
{
    std::puts("[1] Hyperspherical parametrisation");
    std::mt19937 rng(1);
    std::uniform_real_distribution<double> U(0.0, 1.0);

    double worst = 0.0;
    for (std::size_t n = 1; n <= 12; ++n)
        for (int rep = 0; rep < 200; ++rep) {
            std::vector<double> w(n), q(n), back(n);
            double s = 0;
            for (auto& x : w) { x = U(rng) + 1e-3; s += x; }
            for (auto& x : w) x /= s;
            fractionsToHyperspherical(w.data(), n, q.data());
            hypersphericalToFractions(q.data(), n, back.data());
            for (std::size_t k = 0; k < n; ++k) worst = std::max(worst, std::fabs(back[k] - w[k]));
        }
    check(worst < 1e-12, "round trip w -> q -> w, n = 1..12 (max |dw|)", worst);

    // same values as the original GetIntensityParameterNew (Type 3: dimension = P[1], stride 2)
    for (unsigned n = 2; n <= 10; ++n) {
        double maxDiff = 0.0;
        for (int rep = 0; rep < 200; ++rep) {
            std::vector<double> q(n), w(n), P(5 + 2 * n, 0.0);
            q[0] = U(rng);
            for (unsigned k = 1; k + 1 < n; ++k) q[k] = U(rng) * M_PI / 2;
            hypersphericalToFractions(q.data(), n, w.data());
            P[1] = n;
            for (unsigned k = 0; k < n; ++k) P[3 + 2 * k] = q[k];
            for (unsigned k = 0; k < n; ++k)
                maxDiff = std::max(maxDiff, std::fabs(original::GetIntensityParameterNew(P.data(), 3, 3, k + 1) - w[k]));
        }
        char label[96];
        if (n == 6) {
            std::snprintf(label, sizeof label, "n = %u: original differs (bug: last uses phi1, not r)", n);
            std::printf("  %-62s info (%.3g)\n", label, maxDiff);
        } else {
            std::snprintf(label, sizeof label, "n = %u: identical to original GetIntensityParameterNew", n);
            check(maxDiff < 1e-12, label, maxDiff);
        }
    }
}

// ---------------------------------------------------------------------------
// Build both parameter vectors for the same physical model.
struct Case {
    std::vector<ResolutionComponent> res;
    std::vector<LifetimeComponent> life;   // order: p-Ps, absolute, free (as sorted in the original)
    std::vector<double> resQ, lifeQ;       // fraction/intensity parameters
    std::vector<double> sigma, t0;
};

static void buildParameters(const Case& c, double bw, std::vector<double>& pOrig, std::vector<double>& pNew,
                            const LifetimeModel& m)
{
    const int R = c.res.size(), L = c.life.size();
    int nFixed = 0;
    for (const auto& l : c.life) nFixed += (!l.pPs && l.intensityMode != ParMode::Free);
    const double B = 12.0, A = 2.0e6;

    pOrig.assign(5 + 3 * R + 2 * L, 0.0);
    pOrig[0] = B; pOrig[1] = R; pOrig[2] = L; pOrig[3] = A * bw;
    for (int i = 0; i < R; ++i) { pOrig[4 + 3 * i] = c.sigma[i]; pOrig[5 + 3 * i] = c.resQ[i]; pOrig[6 + 3 * i] = c.t0[i]; }
    for (int j = 0; j < L; ++j) { pOrig[4 + 3 * R + 2 * j] = c.life[j].tau; pOrig[5 + 3 * R + 2 * j] = c.lifeQ[j]; }
    pOrig[4 + 3 * R + 2 * L] = nFixed;   // = P[6 + 3R + 2(L-1)]

    pNew.assign(m.nPar(), 0.0);
    pNew[0] = B; pNew[1] = A;
    for (int i = 0; i < R; ++i) { pNew[m.sigmaPar(i)] = c.sigma[i]; pNew[m.offsetPar(i)] = c.t0[i]; pNew[m.resolutionFractionPar(i)] = c.resQ[i]; }
    for (int j = 0; j < L; ++j) { pNew[m.tauPar(j)] = c.life[j].tau; pNew[m.intensityPar(j)] = c.lifeQ[j]; }
}

static LifetimeComponent comp(double tau, ParMode iMode, bool ps = false)
{
    LifetimeComponent c;
    c.tau = tau;
    c.intensity = 0.1;
    c.intensityMode = iMode;
    c.pPs = ps;
    return c;
}

static Case caseNoPs()
{
    Case c;
    c.res = {ResolutionComponent{}, ResolutionComponent{}};
    c.sigma = {0.12, 0.25};
    c.t0 = {6.10, 6.25};
    c.resQ = {0.8, 0.0};                                   // r, unused
    c.life = {comp(0.374, ParMode::Fixed), comp(0.15, ParMode::Free), comp(0.42, ParMode::Free),
              comp(2.5, ParMode::Free), comp(9.0, ParMode::Free)};
    c.lifeQ = {0.08, 0.7, 0.9, 0.4, 0.0};                  // abs I; r, phi1, phi2, unused
    return c;
}

static Case casePs()
{
    Case c;
    c.res = {ResolutionComponent{}};
    c.sigma = {0.15};
    c.t0 = {6.12};
    c.resQ = {0.0};
    c.life = {comp(0.125, ParMode::Fixed, true), comp(1.8, ParMode::Fixed), comp(0.374, ParMode::Fixed),
              comp(0.2, ParMode::Free), comp(0.45, ParMode::Free), comp(3.0, ParMode::Free)};
    c.lifeQ = {1.0 / 3.0, 0.05, 0.1, 0.75, 0.6, 0.0};      // ratio, abs, abs; r, phi1, unused
    return c;
}

static double compareModels(const Case& c, bool ps)
{
    const double bw = 0.02;
    LifetimeModel m(c.res, c.life, bw, 0.7);
    std::vector<double> pO, pN;
    buildParameters(c, bw, pO, pN, m);
    // difference relative to the peak height; the original ExMoGa loses precision
    // in the leading edge (1 + erf(z) for z << 0), hence not exactly 0
    double worst = 0.0, peak = 0.0;
    for (double t = 5.0; t < 40.0; t += 0.01) {
        double x = t;
        const double vo = ps ? original::DiscreteFitFunctionPs(&x, pO.data()) : original::DiscreteFitFunctionNoPs(&x, pO.data());
        const double vn = m(&x, pN.data());
        worst = std::max(worst, std::fabs(vn - vo));
        peak = std::max(peak, std::fabs(vo));
    }
    worst /= peak;
    // fractions sum to 1
    const auto f = m.fractions(pN.data());
    double s = 0.0;
    for (double v : f.lifetime) s += v;
    check(std::fabs(s - 1.0) < 1e-12, ps ? "  p-Ps case: lifetime intensities sum to 1" : "  no-Ps case: lifetime intensities sum to 1", s - 1.0);
    return worst;
}

static void testModel()
{
    std::puts("[2] Full model vs original DiscreteFitFunction*");
    const double d1 = compareModels(caseNoPs(), false);
    check(d1 < 1e-10, "no p-Ps: 2 Gauss, 1 fixed + 4 free (max diff / peak)", d1);
    const double d2 = compareModels(casePs(), true);
    check(d2 < 1e-10, "p-Ps: 1 Gauss, p-Ps + 2 fixed + 3 free (max diff / peak)", d2);
}

// ---------------------------------------------------------------------------
static void testEmg()
{
    std::puts("[3] EMG");
    // direct formula in long double where it is safe
    double worst = 0.0;
    for (double tau : {0.05, 0.125, 0.4, 2.0, 50.0})
        for (double sigma : {0.08, 0.15, 0.3})
            for (double t = 4.0; t < 30.0; t += 0.013) {
                const long double s = sigma, ta = tau, d = t - 6.0;
                const long double arg = s * s / (2 * ta * ta) - d / ta;
                if (arg > 11000) continue;
                const long double ref = std::exp(arg) * std::erfc((s * s / ta - d) / (std::sqrt(2.0L) * s)) / (2 * ta);
                const double v = emg(t, sigma, 6.0, tau);
                if (ref > 1e-300) worst = std::max(worst, (double)std::fabs((v - ref) / ref));
            }
    check(worst < 1e-8, "relative accuracy vs long double formula", worst);

    double worstOrig = 0.0;
    for (double t = 5.0; t < 7.0; t += 0.01) {
        const long double s = 0.12, ta = 0.15, d = t - 6.1;
        const long double ref = std::exp(s * s / (2 * ta * ta) - d / ta) * std::erfc((s * s / ta - d) / (std::sqrt(2.0L) * s)) / (2 * ta);
        worstOrig = std::max(worstOrig, (double)std::fabs((original::ExMoGa(t, 0.12, 6.1, 0.15, 1.0) - ref) / ref));
    }
    std::printf("  %-62s info (%.3g)\n", "original ExMoGa, leading edge: max relative error", worstOrig);

    bool finite = true;
    for (double tau : {0.005, 0.01, 0.05})
        for (double t = -50.0; t < 10.0; t += 0.1) finite = finite && std::isfinite(emg(t, 0.1, 6.0, tau));
    check(finite, "finite for tiny tau / far left tail (original overflows)", 0.0);

    double integral = 0.0;
    for (double t = -10.0; t < 120.0; t += 0.0005) integral += emg(t, 0.15, 6.0, 2.0) * 0.0005;
    check(std::fabs(integral - 1.0) < 1e-6, "unit area", integral - 1.0);
}

// ---------------------------------------------------------------------------
static void benchmark()
{
    std::puts("[4] Benchmark: one fit-like pass = 2000 bins x 400 parameter sets");
    for (bool ps : {false, true}) {
        const Case c = ps ? casePs() : caseNoPs();
        const double bw = 0.02;
        LifetimeModel m(c.res, c.life, bw, 0.7);
        std::vector<double> pO, pN;
        buildParameters(c, bw, pO, pN, m);
        const int nBins = 2000, nSets = 400;
        volatile double sink = 0.0;

        auto t0 = std::chrono::steady_clock::now();
        for (int s = 0; s < nSets; ++s) {
            pO[3] *= 1.0 + 1e-9;
            for (int b = 0; b < nBins; ++b) { double x = 5.0 + 0.02 * b; sink = sink + (ps ? original::DiscreteFitFunctionPs(&x, pO.data()) : original::DiscreteFitFunctionNoPs(&x, pO.data())); }
        }
        auto t1 = std::chrono::steady_clock::now();
        for (int s = 0; s < nSets; ++s) {
            pN[1] *= 1.0 + 1e-9;
            for (int b = 0; b < nBins; ++b) { double x = 5.0 + 0.02 * b; sink = sink + m(&x, pN.data()); }
        }
        auto t2 = std::chrono::steady_clock::now();
        // new model with fractions recomputed in every call (what the cache saves)
        for (int s = 0; s < nSets; ++s)
            for (int b = 0; b < nBins; ++b) { double x = 5.0 + 0.02 * b; pN[1] *= 1.0 + 1e-15; sink = sink + m(&x, pN.data()); }
        auto t3 = std::chrono::steady_clock::now();
        auto ms = [](auto a, auto b) { return std::chrono::duration<double, std::milli>(b - a).count(); };
        std::printf("  %-8s original %7.1f ms | new %7.1f ms (x%.2f) | new without cache %7.1f ms\n",
                    ps ? "p-Ps:" : "no p-Ps:", ms(t0, t1), ms(t1, t2), ms(t0, t1) / ms(t1, t2), ms(t2, t3));
    }

    std::puts("  Fractions only, 10^6 evaluations of all n intensities:");
    std::mt19937 rng(2);
    std::uniform_real_distribution<double> U(0.1, 1.4);
    for (unsigned n : {3u, 5u, 7u}) {
        std::vector<double> q(n), w(n), P(5 + 2 * n, 0.0);
        for (auto& x : q) x = U(rng);
        q[0] = 0.7;
        P[1] = n;
        for (unsigned k = 0; k < n; ++k) P[3 + 2 * k] = q[k];
        volatile double sink = 0.0;
        auto t0 = std::chrono::steady_clock::now();
        for (int r = 0; r < 1000000; ++r) {
            P[3] = 0.7 + 1e-12 * (r & 1);
            for (unsigned k = 1; k <= n; ++k) sink = sink + original::GetIntensityParameterNew(P.data(), 3, 3, k);
        }
        auto t1 = std::chrono::steady_clock::now();
        for (int r = 0; r < 1000000; ++r) {
            q[0] = 0.7 + 1e-12 * (r & 1);
            hypersphericalToFractions(q.data(), n, w.data());
            for (unsigned k = 0; k < n; ++k) sink = sink + w[k];
        }
        auto t2 = std::chrono::steady_clock::now();
        auto ms = [](auto a, auto b) { return std::chrono::duration<double, std::milli>(b - a).count(); };
        std::printf("    n = %u: original switch (per component) %6.1f ms | generic loop (all at once) %6.1f ms\n",
                    n, ms(t0, t1), ms(t1, t2));
    }
}

// ---------------------------------------------------------------------------
// The formula of config/tumor_formula.json, written as a C++ lambda
// (parameter order: A, I1, I3, tau1, tau3, tau4).
static double tumorDecay(double x, const double* q)
{
    const double A = q[0], I1 = q[1], I3 = q[2], t1 = q[3], t3 = q[4], t4 = q[5];
    return A * (I1 * std::exp(-x / t1) / t1 + 0.1 * std::exp(-x / 0.374) / 0.374 + I3 * std::exp(-x / t3) / t3 +
                (0.9 - I1 - I3) * std::exp(-x / t4) / t4);
}

static void testFormulaModel()
{
    std::puts("[5] Formula model vs lifetime model (same physics)");
    const double bw = 0.02, sigma = 0.115, t0 = 6.0, B = 20.0, A = 7.0e5;
    const double I1 = 0.2, I3 = 0.52, t1 = 0.165, t3 = 0.44, t4 = 1.98;

    // lifetime model: 4 absolute intensities (all "fixed" mode -> parameters = intensities)
    std::vector<ResolutionComponent> res(1);
    res[0].sigma = sigma;
    std::vector<LifetimeComponent> life;
    for (double tau : {t1, 0.374, t3, t4}) life.push_back(comp(tau, ParMode::Fixed));
    LifetimeModel lm(res, life, bw, 0.7);
    std::vector<double> pl(lm.nPar(), 0.0);
    pl[0] = B; pl[1] = A; pl[lm.sigmaPar(0)] = sigma; pl[lm.offsetPar(0)] = t0;
    const double I[4] = {I1, 0.1, I3, 0.9 - I1 - I3}, T[4] = {t1, 0.374, t3, t4};
    for (int j = 0; j < 4; ++j) { pl[lm.tauPar(j)] = T[j]; pl[lm.intensityPar(j)] = I[j]; }

    const std::vector<std::string> names = {"A", "I1", "I3", "tau1", "tau3", "tau4"};
    for (int pps : {5, 10}) {
        FormulaModel fm(tumorDecay, names, true, res, bw, sigma / pps, 8.0, 40.0);
        std::vector<double> pf(fm.nPar(), 0.0);
        pf[0] = B; pf[fm.sigmaPar(0)] = sigma; pf[fm.offsetPar(0)] = t0;
        const double q[6] = {A, I1, I3, t1, t3, t4};
        std::copy(q, q + 6, pf.begin() + fm.firstFormulaPar());
        double worst = 0.0, peak = 0.0;
        for (double t = 5.0; t < 45.0; t += 0.02) {   // beyond the initial grid -> grid extension
            const double vl = lm(&t, pl.data()), vf = fm(&t, pf.data());
            worst = std::max(worst, std::fabs(vl - vf));
            peak = std::max(peak, vl - B);
        }
        char label[96];
        std::snprintf(label, sizeof label, "convolve = true, %2d points/sigma (max diff / peak)", pps);
        check(worst / peak < (pps == 5 ? 5e-5 : 5e-6), label, worst / peak);
    }

    // no convolution: formula with pals_emg reproduces the lifetime model exactly
    auto emgFormula = [](double x, const double* q) {
        const double A = q[0], I1 = q[1], I3 = q[2], t1 = q[3], t3 = q[4], t4 = q[5], s = q[6], t0 = q[7];
        return A * (I1 * pals_emg(x, s, t0, t1) + 0.1 * pals_emg(x, s, t0, 0.374) + I3 * pals_emg(x, s, t0, t3) +
                    (0.9 - I1 - I3) * pals_emg(x, s, t0, t4));
    };
    FormulaModel fe(emgFormula, {"A", "I1", "I3", "tau1", "tau3", "tau4", "sigma", "t0"}, false, {}, bw, 0.0, 8.0, 40.0);
    std::vector<double> pe = {B, A, I1, I3, t1, t3, t4, sigma, t0};
    double worst = 0.0;
    for (double t = 5.0; t < 45.0; t += 0.02) worst = std::max(worst, std::fabs(lm(&t, pl.data()) - fe(&t, pe.data())));
    check(worst < 1e-8, "convolve = false with pals_emg (max |diff| in counts)", worst);
}

// ---------------------------------------------------------------------------
// Simulated spectrum: 0.16 ns (25 %), 0.42 ns (55 %), o-Ps 20 % as S0/S1/S2
static LifetimeModel simulationModel(const std::string& scen, std::vector<double>& p, double N)
{
    std::vector<ResolutionComponent> res(1);
    res[0].sigma = 0.115;
    std::vector<LifetimeComponent> life = {comp(0.16, ParMode::Fixed), comp(0.42, ParMode::Fixed)};
    std::vector<double> I = {0.25, 0.55};
    if (scen == "S0") { life.push_back(comp(2.25, ParMode::Fixed)); I.push_back(0.20); }
    if (scen == "S1") { life.push_back(comp(2.0, ParMode::Fixed)); life.push_back(comp(2.5, ParMode::Fixed)); I.push_back(0.1); I.push_back(0.1); }
    if (scen == "S2") { auto c = comp(2.25, ParMode::Fixed); c.hasWidth = true; life.push_back(c); I.push_back(0.20); }
    LifetimeModel m(res, life, 0.02, 0.7);
    p.assign(m.nPar(), 0.0);
    p[0] = 3e-5 * N; p[1] = N; p[m.sigmaPar(0)] = 0.115; p[m.offsetPar(0)] = 6.0;
    for (int j = 0; j < m.nLifetimes(); ++j) { p[m.tauPar(j)] = life[j].tau; p[m.intensityPar(j)] = I[j]; }
    if (scen == "S2") p[m.widthPar(2)] = std::sqrt(std::log(1.0 + std::pow(0.25 / 2.25, 2)));   // sd 0.25 ns
    return m;
}

static void testWidthsAndDistribution()
{
    std::puts("[6] Component width, quadrature, NNLS, distribution analysis");
    GaussHermite gh(32);
    double s0 = 0, s2 = 0, s4 = 0;
    for (std::size_t k = 0; k < gh.x.size(); ++k) { s0 += gh.w[k]; s2 += gh.w[k] * gh.x[k] * gh.x[k]; s4 += gh.w[k] * std::pow(gh.x[k], 4); }
    check(std::fabs(s0 - 1) < 1e-12 && std::fabs(s2 - 1) < 1e-12 && std::fabs(s4 - 3) < 1e-10, "Gauss-Hermite: moments 1, 1, 3 of N(0,1)", s4 - 3);

    // width model vs brute-force integration over a log-normal
    std::vector<double> p;
    LifetimeModel m = simulationModel("S2", p, 1e6);
    const double w = p[m.widthPar(2)], mu = std::log(2.25) - 0.5 * w * w;
    double worst = 0.0, peak = 0.0;
    for (double t = 5.0; t < 40.0; t += 0.05) {
        double ref = 0.0, norm = 0.0;
        for (int k = -4000; k <= 4000; ++k) {
            const double z = k * 0.002, wt = std::exp(-0.5 * z * z);
            ref += wt * emg(t, 0.115, 6.0, std::exp(mu + w * z));
            norm += wt;
        }
        double tt = t;
        const double v = m(&tt, p.data());
        const double other = p[0] + p[1] * 0.02 * (0.25 * emg(t, 0.115, 6.0, 0.16) + 0.55 * emg(t, 0.115, 6.0, 0.42));
        worst = std::max(worst, std::fabs((v - other) - p[1] * 0.02 * 0.2 * ref / norm));
        peak = std::max(peak, v);
    }
    check(worst / peak < 1e-9, "log-normal component = direct integration (max diff / peak)", worst / peak);

    // distribution analysis on simulated spectra, 1e7 counts
    DistributionConfig dc;
    dc.enabled = true;
    std::mt19937_64 rng(11);
    for (const char* scen : {"S0", "S1", "S2"}) {
        LifetimeModel sm = simulationModel(scen, p, 1e7);
        DistributionProblem pr;
        pr.binWidth = 0.02;
        pr.resolution = {{1.0, 0.115, 6.0}};
        pr.background = p[0];
        for (double t = 5.85; t <= 40.0; t += 0.02) {          // fit range like the program
            double tt = t;
            const double mexp = sm(&tt, p.data());
            std::poisson_distribution<long long> pois(mexp);
            pr.t.push_back(t);
            pr.y.push_back(static_cast<double>(pois(rng)));
            pr.weight.push_back(1.0 / std::max(mexp, 1.0));
        }
        auto t0 = std::chrono::steady_clock::now();
        const auto a = analyzeDistribution(pr, dc, {{1.0, 6.0}});
        const double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        const auto& bag = a.estimates[0].regions[0];
        const auto& me = a.estimates[1].regions[0];
        std::printf("  %s: NNLS-bagged I %.2f %% c %.3f w %.3f [%.3f,%.3f] | MaxEnt I %.2f %% c %.3f w %.3f | %.1f s\n", scen,
                    100 * bag.intensity, bag.centroid, bag.width, a.estimates[0].spreads[0].widthQ16,
                    a.estimates[0].spreads[0].widthQ84, 100 * me.intensity, me.centroid, me.width, sec);
        for (const auto& l : a.comparison) std::printf("      %s\n", l.c_str());
        char label[96];
        std::snprintf(label, sizeof label, "%s: o-Ps intensity 20 %% and centroid 2.25 ns recovered (both methods)", scen);
        const bool ok = std::fabs(bag.intensity - 0.2) < 0.01 && std::fabs(me.intensity - 0.2) < 0.01 &&
                        std::fabs(bag.centroid - 2.25) < 0.05 && std::fabs(me.centroid - 2.25) < 0.05;
        check(ok, label, me.centroid);
    }
}

// ---------------------------------------------------------------------------
static void testData(const char* dataPath, const char* configPath)
{
    std::puts("[7] Real spectrum: preliminary estimates and the model at the initial parameters");
    const Config cfg = Config::load(configPath);
    const Spectrum s = binTimes(readTextColumn(dataPath, cfg.input.column), cfg.histogram);
    const Estimates e = estimate(s, cfg.estimation);
    printEstimates(std::cout, s, e);

    auto res = cfg.resolution;
    for (auto& r : res) if (r.sigmaAuto) r.sigma = e.sigma;
    LifetimeModel m(res, cfg.lifetimes, s.binWidth, cfg.fit.oPsThreshold);
    std::vector<double> p;
    m.setInitialIntensityParameters(p);
    p[0] = e.background;
    p[1] = e.area;
    for (int i = 0; i < m.nResolution(); ++i) { p[m.sigmaPar(i)] = res[i].sigma; p[m.offsetPar(i)] = e.peakPosition; }
    for (int j = 0; j < m.nLifetimes(); ++j) p[m.tauPar(j)] = cfg.lifetimes[j].tau;

    double data = 0.0, model = 0.0, chi2 = 0.0;
    int nb = 0;
    for (std::size_t i = e.fitFrom; i <= e.fitTo; ++i) {
        double x = s.center(i);
        const double v = m(&x, p.data());
        data += s.counts[i];
        model += v;
        if (s.counts[i] > 0) { chi2 += (s.counts[i] - v) * (s.counts[i] - v) / s.counts[i]; ++nb; }
    }
    std::printf("  counts in fit range: data %.0f, model(initial) %.0f, ratio %.3f; chi2/bin(initial) %.1f\n",
                data, model, model / data, chi2 / nb);
    check(std::fabs(model / data - 1.0) < 0.2, "area estimate gives the right scale (no 0.013 factor)", model / data);

    // formula configs: parsing and initial parameters
    const Config fc = Config::load("config/tumor_formula.json");
    const auto fres = resolvedResolution(fc, e);
    FormulaModel fm(tumorDecay, {"A", "I1", "I3", "tau1", "tau3", "tau4"}, true, fres, s.binWidth,
                    formulaGridStep(fc, fres), fc.model.windowSigma, 40.0);
    const auto specs = formulaSpecs(fm, fc, s, e);
    check(specs[fm.firstFormulaPar()].value == e.area && specs[fm.firstFormulaPar() + 3].releaseStage == 1,
          "tumor_formula.json: A = estimate:area, tau released in stage 1", specs[fm.firstFormulaPar()].value);
    const Config ec = Config::load("config/tumor_formula_emg.json");
    check(!ec.model.convolve && ec.resolution.empty() && ec.model.parameters.size() == 8,
          "tumor_formula_emg.json: no convolution, 8 parameters", ec.model.parameters.size());

    // config round trip: toJson() of a parsed config parses to the same config
    for (const char* f : {"config/tumor.json", "config/xad4.json", "config/tumor_formula.json",
                          "config/tumor_formula_emg.json", "config/tumor_distribution.json"}) {
        const Config a = Config::load(f);
        const std::string ja = a.toJson();
        const Config b = Config::parse(ja);
        bool same = b.toJson() == ja && a.lifetimes.size() == b.lifetimes.size() && a.resolution.size() == b.resolution.size();
        for (std::size_t j = 0; same && j < a.lifetimes.size(); ++j)
            same = a.lifetimes[j].tau == b.lifetimes[j].tau && a.lifetimes[j].intensity == b.lifetimes[j].intensity &&
                   a.lifetimes[j].tauMode == b.lifetimes[j].tauMode && a.lifetimes[j].intensityMode == b.lifetimes[j].intensityMode &&
                   a.lifetimes[j].pPs == b.lifetimes[j].pPs && a.lifetimes[j].hasWidth == b.lifetimes[j].hasWidth;
        same = same && a.estimation.fitEnd == b.estimation.fitEnd && a.fit.options == b.fit.options &&
               a.model.formula == b.model.formula && a.distribution.enabled == b.distribution.enabled &&
               a.input.format == b.input.format && a.histogram.binWidth == b.histogram.binWidth;
        char label[96];
        std::snprintf(label, sizeof label, "round trip toJson -> parse: %s", f);
        check(same, label, 0.0);
    }

    // explicit fit start and background window (as set with the mouse in the GUI)
    Config r = cfg;
    r.estimation.useFitStart = true;
    r.estimation.fitStart = 5.5;
    r.estimation.fitEnd = 30.0;
    r.estimation.background.useRange = true;
    r.estimation.background.from = 50.0;
    r.estimation.background.to = 100.0;
    const Config r2 = Config::parse(r.toJson());
    const Estimates er = estimate(s, r2.estimation);
    const bool okRange = std::fabs(s.center(er.fitFrom) - 5.51) < 0.011 && std::fabs(s.center(er.fitTo) - 30.01) < 0.011 &&
                         std::fabs(s.center(er.backgroundFrom) - 50.01) < 0.011 && std::fabs(s.center(er.backgroundTo) - 99.99) < 0.011;
    std::printf("  fit %.3f - %.3f ns, background %.3f - %.3f ns: %.3f +- %.3f\n", s.center(er.fitFrom), s.center(er.fitTo),
                s.center(er.backgroundFrom), s.center(er.backgroundTo), er.background, er.backgroundError);
    check(okRange, "explicit fit start / end and background window (after JSON round trip)", er.background);
}

int main(int argc, char** argv)
{
    testHyperspherical();
    testModel();
    testEmg();
    benchmark();
    testFormulaModel();
    testWidthsAndDistribution();
    if (argc == 3) testData(argv[1], argv[2]);
    std::printf("\n%s (%d failure%s)\n", failures ? "FAILED" : "ALL TESTS PASSED", failures, failures == 1 ? "" : "s");
    return failures ? 1 : 0;
}
