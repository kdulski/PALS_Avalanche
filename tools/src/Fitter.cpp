#include "Fitter.h"

#include <TF1.h>
#include <TFitResult.h>
#include <TH1D.h>
#include <TMath.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <stdexcept>

namespace pals {

Fitter::Fitter(const Config& cfg, const Spectrum& spectrum, const Estimates& est, std::vector<ParameterSpec> specs)
    : cfg_(cfg), spectrum_(spectrum), est_(est), specs_(std::move(specs))
{
    const int n = static_cast<int>(spectrum_.size());
    hist_ = std::make_unique<TH1D>("data", "PAL spectrum;time [ns];counts", n, spectrum_.min, spectrum_.max());
    hist_->SetDirectory(nullptr);
    for (int i = 0; i < n; ++i) hist_->SetBinContent(i + 1, spectrum_.counts[i]);
    hist_->SetEntries(std::accumulate(spectrum_.counts.begin(), spectrum_.counts.end(), 0.0));
}

Fitter::~Fitter() = default;

double Fitter::fitMin() const { return spectrum_.center(est_.fitFrom); }
double Fitter::fitMax() const { return spectrum_.center(est_.fitTo); }

void Fitter::initialiseFunction()
{
    func_->SetNpx(std::max<int>(100, static_cast<int>(est_.fitTo - est_.fitFrom + 1)));
    for (std::size_t k = 0; k < specs_.size(); ++k) {
        func_->SetParName(k, specs_[k].name.c_str());
        func_->SetParameter(k, specs_[k].value);
    }
}

// --- parameter state helpers (state is kept in the TF1: fixed <=> limits low >= high)
void Fitter::fix(int ipar)
{
    func_->FixParameter(ipar, func_->GetParameter(ipar));
}

void Fitter::limit(int ipar, double lo, double hi)
{
    func_->ReleaseParameter(ipar);
    func_->SetParameter(ipar, std::clamp(func_->GetParameter(ipar), lo, hi));
    func_->SetParLimits(ipar, lo, hi);
}

void Fitter::partial(int ipar)
{
    const double v = func_->GetParameter(ipar);
    const double d = std::fabs(v) * cfg_.fit.partialVariation;
    if (d > 0.0) limit(ipar, v - d, v + d);
    else fix(ipar);
}

void Fitter::configureStage(int stage)
{
    for (std::size_t k = 0; k < specs_.size(); ++k) {
        const auto& s = specs_[k];
        const int ip = static_cast<int>(k);
        if (s.mode == ParMode::Fixed || stage < s.releaseStage || (s.mode == ParMode::Partial && stage == 0))
            fix(ip);
        else if (s.mode == ParMode::Partial)
            partial(ip);
        else if (s.hasLimits)
            limit(ip, s.lower, s.upper);
        else
            func_->ReleaseParameter(ip);
    }
}

StageSummary Fitter::fitStage(int stage, const std::string& options)
{
    configureStage(stage);
    std::cout << "\n------------- Fit stage " << stage << " (options \"" << options << "\") -------------\n";
    // S: return the fit result, 0: do not draw (the function is still stored in the histogram)
    result_ = hist_->Fit(func_.get(), (options + "S0").c_str());
    if (!result_.Get()) throw std::runtime_error("Fit stage " + std::to_string(stage) + " failed (no fit result)");

    StageSummary s;
    s.stage = stage;
    s.options = options;
    s.status = result_->Status();
    s.chi2 = func_->GetChisquare();
    s.ndf = func_->GetNDF();
    return s;
}

FitOutcome Fitter::run()
{
    stages_.clear();
    stages_.push_back(fitStage(0, cfg_.fit.firstPassOptions));
    if (onStage_) onStage_(stages_.back());
    for (int it = 1; it <= cfg_.fit.iterations; ++it) {
        stages_.push_back(fitStage(it, cfg_.fit.options));
        if (onStage_) onStage_(stages_.back());
    }
    FitOutcome out = collect();
    out.objective = objectiveName();

    std::vector<int> pars;
    for (int k : profilePars_)
        if (!out.fixed[k]) pars.push_back(k);
    if (pars.empty()) return out;
    if (stages_.back().options.find('W') != std::string::npos) {
        std::cerr << "Note: no profile intervals with option W (unit weights)\n";
        return out;
    }
    for (int k : pars) {
        std::cout << "\n------------- Profile likelihood of " << out.names[k] << " -------------\n";
        out.profiles.push_back(profile(k, out));
    }
    // restore the best fit (function stored in the histogram, fit result)
    applyState(out);
    func_->SetParameters(out.params.data());
    result_ = hist_->Fit(func_.get(), (stages_.back().options + "QS0").c_str());
    return out;
}

// Statistic minimised by ROOT for the options of the last stage, recomputed here
// (Neyman chi2 with sqrt(N) errors and empty bins skipped; P: Pearson chi2;
//  L: Poisson deviance = Baker-Cousins likelihood chi2). In all cases a rise by 1
// gives the 68 % interval of one parameter.
std::string Fitter::objectiveName() const
{
    const std::string& o = stages_.back().options;
    if (o.find('L') != std::string::npos) return "Poisson deviance";
    if (o.find('P') != std::string::npos) return "Pearson chi2";
    return "Neyman chi2";
}

double Fitter::objective() const
{
    const std::string& o = stages_.back().options;
    const bool poisson = o.find('L') != std::string::npos;
    const bool pearson = o.find('P') != std::string::npos;
    double sum = 0.0;
    for (std::size_t i = est_.fitFrom; i <= est_.fitTo; ++i) {
        const double y = spectrum_.counts[i];
        const double m = std::max(func_->Eval(spectrum_.center(i)), 1e-12);
        if (poisson) sum += 2.0 * (m - y + (y > 0.0 ? y * std::log(y / m) : 0.0));
        else if (pearson) sum += (y - m) * (y - m) / m;
        else if (y > 0.0) sum += (y - m) * (y - m) / y;
    }
    return sum;
}

// fixed / limited / free exactly as after the last stage (no re-centring of windows)
void Fitter::applyState(const FitOutcome& state)
{
    for (std::size_t k = 0; k < specs_.size(); ++k) {
        const int ip = static_cast<int>(k);
        if (state.fixed[k]) fix(ip);
        else if (state.lower[k] == 0.0 && state.upper[k] == 0.0) func_->ReleaseParameter(ip);
        else limit(ip, state.lower[k], state.upper[k]);
    }
}

ProfileInterval Fitter::profile(int ipar, const FitOutcome& base)
{
    ProfileInterval pi;
    pi.par = ipar;
    pi.best = base.params[ipar];
    const std::string opt = stages_.back().options + "QS0";
    func_->SetParameters(base.params.data());
    const double f0 = objective();
    const bool bounded = !(base.lower[ipar] == 0.0 && base.upper[ipar] == 0.0);
    const double lo = bounded ? base.lower[ipar] : -1e300, hi = bounded ? base.upper[ipar] : 1e300;
    const double step = 0.5 * (base.errors[ipar] > 0.0 ? base.errors[ipar] : 0.02 * std::max(1.0, std::fabs(pi.best)));

    for (int dir : {-1, +1}) {
        func_->SetParameters(base.params.data());
        double prevV = pi.best, prevD = 0.0;
        double bound = dir < 0 ? lo : hi;
        bool atLimit = true;
        for (int k = 1; k <= 40; ++k) {
            double v = pi.best + dir * step * k;
            const bool clipped = dir < 0 ? v <= lo : v >= hi;
            if (clipped) v = dir < 0 ? lo : hi;
            applyState(base);
            func_->FixParameter(ipar, v);
            hist_->Fit(func_.get(), opt.c_str());
            ++pi.fits;
            const double d = objective() - f0;
            if (d < -0.05) pi.lowerMinimumFound = true;
            if (d >= 1.0) {
                bound = prevV + (v - prevV) * (1.0 - prevD) / (d - prevD);
                atLimit = false;
                break;
            }
            prevV = v;
            prevD = d;
            if (clipped) break;
        }
        if (atLimit) bound = prevV;      // limit reached (or 20 sigma) without a rise by 1
        if (dir < 0) { pi.lower = bound; pi.lowerAtLimit = atLimit; }
        else { pi.upper = bound; pi.upperAtLimit = atLimit; }
    }
    return pi;
}

FitOutcome Fitter::collect() const
{
    const int n = static_cast<int>(specs_.size());
    FitOutcome out;
    out.params.resize(n);
    out.errors.resize(n);
    out.lower.resize(n);
    out.upper.resize(n);
    out.fixed.resize(n);
    out.covariance.assign(n, std::vector<double>(n, 0.0));
    for (int k = 0; k < n; ++k) {
        out.names.push_back(specs_[k].name);
        out.params[k] = func_->GetParameter(k);
        double lo = 0.0, hi = 0.0;
        func_->GetParLimits(k, lo, hi);
        out.lower[k] = lo;
        out.upper[k] = hi;
        out.fixed[k] = result_->IsParameterFixed(k);
        out.errors[k] = out.fixed[k] ? 0.0 : func_->GetParError(k);
    }
    for (int a = 0; a < n; ++a)
        for (int b = 0; b < n; ++b) out.covariance[a][b] = result_->CovMatrix(a, b);

    out.chi2 = func_->GetChisquare();
    out.ndf = func_->GetNDF();
    out.pValue = out.ndf > 0 ? TMath::Prob(out.chi2, out.ndf) : 0.0;
    out.status = result_->Status();
    out.covarianceStatus = result_->CovMatrixStatus();
    out.fitMin = fitMin();
    out.fitMax = fitMax();
    out.fitBins = static_cast<int>(est_.fitTo - est_.fitFrom + 1);
    out.stages = stages_;
    return out;
}

} // namespace pals
