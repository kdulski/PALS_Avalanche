#include "TextReport.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <numeric>
#include <ostream>
#include <sstream>

namespace pals {

std::vector<ValueError> propagate(const std::function<std::vector<double>(const double*)>& g, const FitOutcome& fit)
{
    const std::size_t nPar = fit.params.size();
    const std::vector<double> central = g(fit.params.data());
    const std::size_t nOut = central.size();

    std::vector<std::vector<double>> jac(nOut, std::vector<double>(nPar, 0.0));
    std::vector<double> p = fit.params;
    for (std::size_t k = 0; k < nPar; ++k) {
        if (fit.fixed[k] || fit.errors[k] <= 0.0) continue;
        const double h = 1e-4 * fit.errors[k];
        p[k] = fit.params[k] + h;
        const auto up = g(p.data());
        p[k] = fit.params[k] - h;
        const auto down = g(p.data());
        p[k] = fit.params[k];
        for (std::size_t o = 0; o < nOut; ++o) jac[o][k] = (up[o] - down[o]) / (2.0 * h);
    }

    std::vector<ValueError> out(nOut);
    for (std::size_t o = 0; o < nOut; ++o) {
        double var = 0.0;
        for (std::size_t a = 0; a < nPar; ++a) {
            if (jac[o][a] == 0.0) continue;
            for (std::size_t b = 0; b < nPar; ++b) var += jac[o][a] * fit.covariance[a][b] * jac[o][b];
        }
        out[o] = {central[o], std::sqrt(std::max(0.0, var))};
    }
    return out;
}

DerivedFractions derivedFractions(const LifetimeModel& model, const FitOutcome& fit)
{
    const auto all = propagate([&](const double* p) {
        const auto f = model.fractions(p);
        std::vector<double> v = f.resolution;
        v.insert(v.end(), f.lifetime.begin(), f.lifetime.end());
        return v;
    }, fit);
    DerivedFractions d;
    d.resolution.assign(all.begin(), all.begin() + model.nResolution());
    d.lifetime.assign(all.begin() + model.nResolution(), all.end());
    return d;
}

namespace {

std::string num(double v, int prec = 5)
{
    std::ostringstream o;
    o << std::setprecision(prec) << v;
    return o.str();
}

std::string withError(double v, double e, int prec = 5) { return num(v, prec) + " +- " + num(e, 2); }

std::string modeOf(const LifetimeComponent& c)
{
    if (c.pPs) return std::string("p-Ps, tau ") + toString(c.tauMode);
    return std::string("tau ") + toString(c.tauMode) + ", I " + toString(c.intensityMode);
}

// TMinuit status = migrad + 10*minos + 100*hesse + 1000*improve.
// MIGRAD 4 = call limit reached / abnormal termination;
// IMPROVE 4 = IMPROVE (option M) found no better minimum - the fit is fine.
std::string statusMeaning(int status)
{
    if (status < 0) return "fit failed";
    const int migrad = status % 10, minos = (status / 10) % 10, hesse = (status / 100) % 10,
              improve = (status / 1000) % 10;
    std::vector<std::string> parts;
    if (migrad == 4) parts.push_back("MIGRAD reached the call limit or terminated abnormally");
    else if (migrad != 0) parts.push_back("MIGRAD did not converge (code " + std::to_string(migrad) + ")");
    else parts.push_back("MIGRAD converged");
    if (minos != 0) parts.push_back("MINOS problem (code " + std::to_string(minos) + ")");
    if (hesse != 0) parts.push_back("HESSE failed, errors unreliable (code " + std::to_string(hesse) + ")");
    if (improve == 4) parts.push_back("IMPROVE found no better minimum");
    else if (improve != 0) parts.push_back("IMPROVE code " + std::to_string(improve));
    std::string out;
    for (std::size_t i = 0; i < parts.size(); ++i) out += (i ? "; " : "") + parts[i];
    return out;
}

std::string covarianceMeaning(int status)
{
    switch (status) {
        case 0: return "not calculated";
        case 1: return "approximate";
        case 2: return "forced positive definite";
        case 3: return "accurate";
        default: return "unknown";
    }
}

void writeResolutionTable(std::ostream& os, const ResolutionLayout& layout, const std::vector<ValueError>& weights,
                          const FitOutcome& fit)
{
    const double fwhm = 2.0 * std::sqrt(2.0 * std::log(2.0));
    const auto& p = fit.params;
    const auto& e = fit.errors;
    os << "Resolution components\n"
       << std::left << std::setw(4) << "#" << std::setw(24) << "sigma [ns]" << std::setw(24) << "FWHM [ns]"
       << std::setw(24) << "t0 [ns]" << std::setw(26) << "fraction" << "modes\n";
    for (int i = 0; i < layout.size(); ++i) {
        const auto& c = layout.components[i];
        const int s = layout.sigmaPar[i], o = layout.offsetPar[i];
        os << std::setw(4) << i + 1 << std::setw(24) << withError(p[s], e[s])
           << std::setw(24) << withError(fwhm * p[s], fwhm * e[s]) << std::setw(24) << withError(p[o], e[o], 6)
           << std::setw(26) << withError(weights[i].value, weights[i].error)
           << "sigma " << toString(c.sigmaMode) << ", fraction " << toString(c.intensityMode) << "\n";
    }
    os << std::right;
}

} // namespace

std::string statisticsNote(double counts, double low, double good, const std::string& what)
{
    std::ostringstream o;
    o << "NOTE on statistics (" << what << "): in simulations of PAL spectra (resolution FWHM 0.27 ns, o-Ps 20 %:\n"
      << "  one component 2.25 ns / 2.0 + 2.5 ns / log-normal with the same mean and width) a broadened o-Ps was\n"
      << "  distinguished from a single exponential reliably only from about " << num(low, 2) << " - " << num(good, 2)
      << " counts.\n"
      << "  Below that, widths are mostly set by noise and regularisation. Two close components and a continuous\n"
      << "  distribution with the same mean and width cannot be told apart (~1e10 counts would be needed).\n"
      << "  This spectrum: " << num(counts, 3) << " counts -> ";
    if (counts < low) o << "BELOW the threshold: treat widths and shapes as not determined by the data.\n";
    else if (counts < good) o << "in the transition range: widths are indicative, check with simulations.\n";
    else o << "above the threshold.\n";
    return o.str();
}

void writeReport(std::ostream& os, const ReportInfo& info, const Spectrum& spectrum, const Estimates& est,
                 const FitOutcome& fit, const ResultsSection& results)
{
    os << "PALS_avalanche - fit of a positron annihilation lifetime spectrum\n"
       << "Date:   " << info.date << "\n"
       << "Data:   " << info.dataPath << "\n"
       << "Config: " << info.configPath << "\n\n";

    os << "=== Preliminary estimates ===\n";
    printEstimates(os, spectrum, est);

    os << "\n=== Fit ===\n"
       << "Range: " << num(fit.fitMin) << " .. " << num(fit.fitMax) << " ns (" << fit.fitBins << " bins)\n";
    for (const auto& s : fit.stages) {
        os << "  stage " << s.stage << "  options " << std::left << std::setw(6) << s.options << std::right
           << "  status " << std::setw(4) << s.status << "  chi2/ndf " << num(s.chi2, 7) << "/" << s.ndf << " = "
           << num(s.ndf > 0 ? s.chi2 / s.ndf : 0.0);
        if (s.options.find('W') != std::string::npos) os << "  (W: unit weights, not a chi2)";
        os << "\n";
    }
    os << "Chi2 / NDF:        " << num(fit.chi2, 7) << " / " << fit.ndf << " = "
       << num(fit.ndf > 0 ? fit.chi2 / fit.ndf : 0.0) << "\n"
       << "p-value:           " << num(fit.pValue, 4) << "\n"
       << "  (probability of a chi2 at least this large from statistical fluctuations alone if the model\n"
       << "   were correct and the bin errors Poissonian; it is NOT the probability that the model is correct.\n"
       << "   Very small values with chi2/NDF close to 1 are common for spectra with ~1e6 counts or more:\n"
       << "   small model or resolution deficiencies become significant.)\n"
       << "Minimizer status:  " << fit.status << " (" << statusMeaning(fit.status) << ")\n"
       << "Covariance matrix: status " << fit.covarianceStatus << " (" << fit.covarianceStatus << " = "
       << covarianceMeaning(fit.covarianceStatus) << ")\n";

    // parameters at their limits
    for (std::size_t k = 0; k < fit.params.size(); ++k) {
        if (fit.fixed[k] || (fit.lower[k] == 0.0 && fit.upper[k] == 0.0)) continue;
        const double tol = 1e-4 * std::max(1e-12, fit.upper[k] - fit.lower[k]);
        if (fit.params[k] - fit.lower[k] < tol || fit.upper[k] - fit.params[k] < tol)
            os << "NOTE: " << fit.names[k] << " is at its limit [" << num(fit.lower[k]) << ", " << num(fit.upper[k]) << "]\n";
    }

    os << "\n=== Results ===\n";
    results(os);

    os << "\n=== All fit parameters ===\n";
    for (std::size_t k = 0; k < fit.params.size(); ++k) {
        os << std::left << std::setw(4) << k << std::setw(26) << fit.names[k] << std::right
           << std::setw(14) << num(fit.params[k], 7) << std::setw(14) << num(fit.errors[k], 3) << "   ";
        if (fit.fixed[k]) os << "fixed";
        else if (fit.lower[k] == 0.0 && fit.upper[k] == 0.0) os << "no limits";
        else os << "[" << num(fit.lower[k]) << ", " << num(fit.upper[k]) << "]";
        os << "\n";
    }
}

void writeLifetimeResults(std::ostream& os, const LifetimeModel& model, const Estimates& est, const FitOutcome& fit)
{
    const auto derived = derivedFractions(model, fit);
    const auto& p = fit.params;
    const auto& e = fit.errors;

    os << "Background [counts/bin]  " << withError(p[LifetimeModel::kBackground], e[LifetimeModel::kBackground], 7)
       << "   (estimate " << withError(est.background, est.backgroundError, 7) << ")\n"
       << "Area [counts]            " << withError(p[LifetimeModel::kArea], e[LifetimeModel::kArea], 8) << "\n\n";

    writeResolutionTable(os, model.resolutionLayout(), derived.resolution, fit);

    std::vector<int> order(model.nLifetimes());
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](int a, int b) { return p[model.tauPar(a)] < p[model.tauPar(b)]; });
    bool anyWidth = false;
    for (int j = 0; j < model.nLifetimes(); ++j) anyWidth = anyWidth || model.widthPar(j) >= 0;

    // intensities renormalised to the components that are not fully fixed
    const auto& comps = model.lifetimes();
    auto fullyFixed = [&](int j) {
        const auto& c = comps[j];
        return !c.pPs && c.tauMode == ParMode::Fixed && c.intensityMode == ParMode::Fixed &&
               (!c.hasWidth || c.widthMode == ParMode::Fixed);
    };
    const auto withoutFixed = propagate([&](const double* q) {
        const auto f = model.fractions(q);
        double rest = 0.0;
        for (int j = 0; j < model.nLifetimes(); ++j) if (!fullyFixed(j)) rest += f.lifetime[j];
        std::vector<double> v(model.nLifetimes(), 0.0);
        for (int j = 0; j < model.nLifetimes(); ++j)
            if (!fullyFixed(j) && rest > 0.0) v[j] = f.lifetime[j] / rest;
        return v;
    }, fit);

    os << "\nLifetime components (sorted by lifetime; # = position in the config)\n"
       << std::left << std::setw(4) << "#" << std::setw(26) << "tau [ns]" << std::setw(26) << "intensity [%]"
       << std::setw(30) << "intensity without fixed [%]";
    if (anyWidth) os << std::setw(42) << "width = sd(ln tau)  [68 % profile]";
    os << "modes\n";
    double total = 0.0;
    for (int j : order) {
        total += derived.lifetime[j].value;
        os << std::setw(4) << j + 1 << std::setw(26) << withError(p[model.tauPar(j)], e[model.tauPar(j)])
           << std::setw(26) << withError(100.0 * derived.lifetime[j].value, 100.0 * derived.lifetime[j].error)
           << std::setw(30)
           << (fullyFixed(j) ? std::string("-")
                             : withError(100.0 * withoutFixed[j].value, 100.0 * withoutFixed[j].error));
        if (anyWidth) {
            std::string w = "-";
            const int iw = model.widthPar(j);
            if (iw >= 0) {
                w = withError(p[iw], e[iw], 3);
                for (const auto& pi : fit.profiles)
                    if (pi.par == iw)
                        w += "  [" + (pi.lowerAtLimit ? std::string("<") : std::string("")) + num(pi.lower, 3) + ", " +
                             num(pi.upper, 3) + (pi.upperAtLimit ? std::string("+") : std::string("")) + "]";
                if (fit.fixed[iw]) w = num(p[iw], 3) + " (fixed)";
            }
            os << std::setw(42) << w;
        }
        os << modeOf(model.lifetimes()[j]) << "\n";
    }
    os << std::right << "Sum of intensities: " << num(100.0 * total) << " %\n";
    if (model.pPsIndex() >= 0)
        os << "p-Ps/o-Ps ratio " << num(p[model.intensityPar(model.pPsIndex())])
           << ", o-Ps = components with tau > " << model.oPsThreshold() << " ns\n";
    if (anyWidth) {
        os << "Width: log-normal distribution of tau with mean tau; sd(tau) ~ tau * sqrt(exp(width^2) - 1).\n";
        if (!fit.profiles.empty())
            os << "Profile intervals: " << fit.objective << " + 1; '<' / '+' = the limit of the parameter was "
               << "reached (the interval extends to it).\n";
        for (const auto& pi : fit.profiles)
            if (pi.lowerMinimumFound) os << "NOTE: a lower minimum was met while profiling " << fit.names[pi.par] << "\n";
        os << statisticsNote(p[LifetimeModel::kArea], 5e6, 1e7, "the width of a component") ;
    }
    os << "'intensity without fixed': intensities renormalised to 100 % over the components whose lifetime\n"
       << "  and intensity (and width) are not all fixed; '-' marks the fully fixed ones.\n"
       << "Intensity errors are propagated with the full covariance matrix.\n";
}

void writeFormulaResults(std::ostream& os, const FormulaModel& model, const std::string& formula,
                         const Estimates& est, const FitOutcome& fit)
{
    const auto& p = fit.params;
    const auto& e = fit.errors;
    os << "Formula: " << formula << "\n";
    if (model.convolve())
        os << "Convolved numerically with the resolution function (grid step " << num(model.step(), 4) << " ns)\n";
    else
        os << "No convolution with the resolution function (it has to be part of the formula)\n";
    os << "\nBackground [counts/bin]  " << withError(p[FormulaModel::kBackground], e[FormulaModel::kBackground], 7)
       << "   (estimate " << withError(est.background, est.backgroundError, 7) << ")\n\n";

    if (model.convolve()) {
        const auto& layout = model.resolutionLayout();
        const auto weights = propagate([&](const double* q) { return layout.weights(q); }, fit);
        writeResolutionTable(os, layout, weights, fit);
        os << "\n";
    }
    os << "Formula parameters\n";
    for (int k = 0; k < model.nFormulaPars(); ++k) {
        const int i = model.firstFormulaPar() + k;
        os << "  " << std::left << std::setw(20) << model.formulaParNames()[k] << std::right
           << withError(p[i], e[i], 7) << (fit.fixed[i] ? "   (fixed)" : "") << "\n";
    }
}

} // namespace pals
