#include "ModelSetup.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace pals {

std::vector<ResolutionComponent> resolvedResolution(const Config& cfg, const Estimates& est)
{
    auto res = cfg.resolution;
    for (auto& r : res)
        if (r.sigmaAuto) r.sigma = est.sigma;
    return res;
}

double estimateValue(const std::string& name, const Spectrum& s, const Estimates& est)
{
    if (name == "background") return est.background;
    if (name == "background_error") return est.backgroundError;
    if (name == "area") return est.area;
    if (name == "peak_position") return est.peakPosition;
    if (name == "peak_counts") return est.peakCounts;
    if (name == "sigma") return est.sigma;
    if (name == "bin_width") return s.binWidth;
    if (name == "fit_start") return s.center(est.fitFrom);
    if (name == "fit_end") return s.center(est.fitTo);
    throw std::runtime_error("Unknown estimate \"" + name + "\" (background, background_error, area, "
                             "peak_position, peak_counts, sigma, bin_width, fit_start, fit_end)");
}

namespace {

ParameterSpec backgroundSpec(const Config& cfg, const Estimates& est)
{
    const double window = cfg.estimation.background.nSigma * est.backgroundError;
    if (window > 0.0) return makeSpec("background", est.background, est.background - window, est.background + window);
    return makeSpec("background", est.background, ParMode::Fixed);
}

// sigma, t0 and fraction parameters of the resolution function
void addResolutionSpecs(std::vector<ParameterSpec>& specs, const ResolutionLayout& layout, const Config& cfg,
                        const Estimates& est, const std::vector<std::string>& names)
{
    std::vector<double> p(specs.size(), 0.0);
    layout.initialFractionParameters(p);
    for (int i = 0; i < layout.size(); ++i) {
        const auto& c = layout.components[i];
        const double sigma = std::clamp(c.sigma, cfg.fit.sigmaMin, cfg.fit.sigmaMax);
        specs[layout.sigmaPar[i]] = makeSpec(names[layout.sigmaPar[i]], sigma, cfg.fit.sigmaMin, cfg.fit.sigmaMax, c.sigmaMode);
        const double w = cfg.fit.offsetWindow * sigma;
        specs[layout.offsetPar[i]] = makeSpec(names[layout.offsetPar[i]], est.peakPosition,
                                              est.peakPosition - w, est.peakPosition + w);
        if (c.intensityMode != ParMode::Free)
            specs[layout.fractionPar[i]] = makeSpec(names[layout.fractionPar[i]], p[layout.fractionPar[i]], 0.0, 1.0,
                                                    c.intensityMode);
    }
    // hyperspherical group: r in [0,1], angles in [0, pi/2], last one unused
    const auto group = layout.groupPars();
    for (std::size_t k = 0; k < group.size(); ++k) {
        const int ip = group[k];
        if (k + 1 == group.size()) specs[ip] = makeSpec(names[ip], p[ip], ParMode::Fixed);
        else specs[ip] = makeSpec(names[ip], p[ip], 0.0, k == 0 ? 1.0 : 0.5 * M_PI);
    }
}

void clampInitial(ParameterSpec& s)
{
    if (!s.hasLimits || (s.value >= s.lower && s.value <= s.upper)) return;
    const double v = std::clamp(s.value, s.lower, s.upper);
    std::cerr << "Warning: initial value of " << s.name << " (" << s.value << ") moved into its limits: " << v << "\n";
    s.value = v;
}

} // namespace

std::vector<ParameterSpec> lifetimeSpecs(const LifetimeModel& model, const Config& cfg, const Estimates& est)
{
    const auto names = model.parameterNames();
    std::vector<ParameterSpec> specs(model.nPar());
    specs[LifetimeModel::kBackground] = backgroundSpec(cfg, est);
    specs[LifetimeModel::kArea] = makeSpec(names[LifetimeModel::kArea], est.area);
    addResolutionSpecs(specs, model.resolutionLayout(), cfg, est, names);

    std::vector<double> p;
    model.setInitialIntensityParameters(p);
    // lifetimes are fixed in stage 0 and released from stage 1
    for (int j = 0; j < model.nLifetimes(); ++j) {
        const auto& c = model.lifetimes()[j];
        const int it = model.tauPar(j), ii = model.intensityPar(j);
        specs[it] = c.tauMode == ParMode::Free
                        ? makeSpec(names[it], c.tau, cfg.fit.tauMin, cfg.fit.tauMax, ParMode::Free, 1)
                        : makeSpec(names[it], c.tau, c.tauMode, 1);
        clampInitial(specs[it]);
        if (c.pPs) specs[ii] = makeSpec(names[ii], p[ii], ParMode::Fixed);
        else if (c.intensityMode != ParMode::Free) specs[ii] = makeSpec(names[ii], p[ii], 0.0, 1.0, c.intensityMode, 1);

        // width (sd of ln tau): released together with the lifetimes; at exactly 0 the
        // model does not depend on it to first order, so a free width starts at >= 0.02
        const int iw = model.widthPar(j);
        if (iw >= 0) {
            double w = c.width;
            if (c.widthMode != ParMode::Fixed && w < 0.02) {
                std::cerr << "Note: initial width of component " << j + 1 << " set to 0.05 (a free width cannot start at 0)\n";
                w = 0.05;
            }
            specs[iw] = makeSpec(names[iw], w, cfg.fit.widthMin, cfg.fit.widthMax, c.widthMode, 1);
            clampInitial(specs[iw]);
        }
    }
    const auto group = model.lifetimeGroupPars();
    for (std::size_t k = 0; k < group.size(); ++k) {
        const int ip = group[k];
        if (k + 1 == group.size()) specs[ip] = makeSpec(names[ip], p[ip], ParMode::Fixed);
        else specs[ip] = makeSpec(names[ip], p[ip], 0.0, k == 0 ? 1.0 : 0.5 * M_PI);
    }
    return specs;
}

double formulaGridStep(const Config& cfg, const std::vector<ResolutionComponent>& resolution)
{
    if (cfg.model.step > 0.0) return cfg.model.step;
    double sigma = 0.0;
    for (const auto& r : resolution) sigma = sigma == 0.0 ? r.sigma : std::min(sigma, r.sigma);
    if (sigma <= 0.0) throw std::runtime_error("Cannot set the convolution grid step (no resolution sigma)");
    return sigma / cfg.model.pointsPerSigma;
}

std::vector<ParameterSpec> formulaSpecs(const FormulaModel& model, const Config& cfg, const Spectrum& s,
                                        const Estimates& est)
{
    const auto names = model.parameterNames();
    std::vector<ParameterSpec> specs(model.nPar());
    specs[FormulaModel::kBackground] = backgroundSpec(cfg, est);
    if (model.convolve()) addResolutionSpecs(specs, model.resolutionLayout(), cfg, est, names);

    const auto& formulaNames = model.formulaParNames();
    for (const auto& par : cfg.model.parameters)
        if (std::find(formulaNames.begin(), formulaNames.end(), par.name) == formulaNames.end())
            std::cerr << "Config warning: parameter \"" << par.name << "\" is not used in the formula\n";

    for (int k = 0; k < model.nFormulaPars(); ++k) {
        const std::string& name = formulaNames[k];
        auto it = std::find_if(cfg.model.parameters.begin(), cfg.model.parameters.end(),
                               [&](const FormulaParameter& f) { return f.name == name; });
        if (it == cfg.model.parameters.end())
            throw std::runtime_error("Formula parameter [" + name + "] is not defined in model.parameters");
        const double value = it->estimate.empty() ? it->value : it->estimateFactor * estimateValue(it->estimate, s, est);
        ParameterSpec spec = it->hasLimits ? makeSpec(name, value, it->lower, it->upper, it->mode, it->stage)
                                           : makeSpec(name, value, it->mode, it->stage);
        clampInitial(spec);
        if (spec.releaseStage > cfg.fit.iterations && spec.mode != ParMode::Fixed)
            std::cerr << "Warning: parameter " << name << " has stage " << spec.releaseStage
                      << " > fit.iterations and stays fixed\n";
        specs[model.firstFormulaPar() + k] = spec;
    }
    return specs;
}

DistributionProblem distributionProblem(const LifetimeModel& model, const FitOutcome& fit, const Spectrum& s,
                                        const Estimates& est)
{
    DistributionProblem pr;
    pr.binWidth = s.binWidth;
    const double* p = fit.params.data();
    const auto w = model.fractions(p).resolution;
    for (int i = 0; i < model.nResolution(); ++i)
        pr.resolution.push_back({w[i], p[model.sigmaPar(i)], p[model.offsetPar(i)]});
    pr.background = p[LifetimeModel::kBackground];
    for (std::size_t i = est.fitFrom; i <= est.fitTo; ++i) {
        double t = s.center(i);
        pr.t.push_back(t);
        pr.y.push_back(s.counts[i]);
        pr.weight.push_back(1.0 / std::max(model(&t, p), 1.0));
    }
    return pr;
}

std::vector<std::pair<double, double>> distributionRegions(const Config& cfg)
{
    if (!cfg.distribution.regions.empty()) return cfg.distribution.regions;
    return {{cfg.fit.oPsThreshold, cfg.distribution.gridMax}};
}

} // namespace pals
