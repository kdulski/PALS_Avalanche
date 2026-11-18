#include "LifetimeModel.h"
#include "Emg.h"
#include "Intensities.h"
#include "Quadrature.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace pals {

double emg(double t, double sigma, double t0, double tau)
{
    return pals_emg(t, sigma, t0, tau);
}

LifetimeModel::LifetimeModel(std::vector<ResolutionComponent> resolution,
                             std::vector<LifetimeComponent> lifetimes,
                             double binWidth, double oPsThreshold)
    : resolution_(std::move(resolution)), lifetimes_(std::move(lifetimes)),
      binWidth_(binWidth), oPsThreshold_(oPsThreshold)
{
    if (resolution_.empty() || lifetimes_.empty())
        throw std::invalid_argument("LifetimeModel: at least one resolution and one lifetime component is needed");

    resLayout_.components = resolution_;
    for (int i = 0; i < nResolution(); ++i) {
        resLayout_.sigmaPar.push_back(sigmaPar(i));
        resLayout_.offsetPar.push_back(offsetPar(i));
        resLayout_.fractionPar.push_back(resolutionFractionPar(i));
    }

    widthPar_.assign(nLifetimes(), -1);
    for (int j = 0; j < nLifetimes(); ++j) {
        if (lifetimes_[j].hasWidth) {
            if (lifetimes_[j].pPs) throw std::invalid_argument("LifetimeModel: the p-Ps component cannot have a width");
            widthPar_[j] = 2 + 3 * nResolution() + 2 * nLifetimes() + nWidths_++;
        }
    }
    for (int j = 0; j < nLifetimes(); ++j) {
        if (lifetimes_[j].pPs) {
            if (pPs_ >= 0) throw std::invalid_argument("LifetimeModel: only one p-Ps component is allowed");
            pPs_ = j;
        } else {
            (lifetimes_[j].intensityMode == ParMode::Free ? lifeFree_ : lifeAbs_).push_back(j);
        }
    }
}

std::vector<int> LifetimeModel::lifetimeGroupPars() const
{
    std::vector<int> pars;
    for (int j : lifeFree_) pars.push_back(intensityPar(j));
    return pars;
}

void LifetimeModel::computeFractions(const double* p, Fractions& f) const
{
    std::vector<double> q, w;

    // --- resolution: absolute ones first, the free group shares the rest
    f.resolution = resLayout_.weights(p);

    // --- lifetimes
    // Every o-Ps component (tau > oPsThreshold) with intensity I is accompanied
    // by p-Ps with intensity ratio*I, where ratio = parameter of the p-Ps component.
    f.lifetime.assign(nLifetimes(), 0.0);
    const bool withPs = pPs_ >= 0;
    const double ratio = withPs ? p[intensityPar(pPs_)] : 0.0;
    auto isOPs = [&](int j) { return withPs && p[tauPar(j)] > oPsThreshold_; };

    double absSum = 0.0;
    double pPsFromAbs = 0.0;
    for (int j : lifeAbs_) {
        const double I = p[intensityPar(j)];
        f.lifetime[j] = I;
        absSum += I;
        if (isOPs(j)) pPsFromAbs += ratio * I;
    }
    absSum += pPsFromAbs;

    double pPsFromFree = 0.0;
    if (!lifeFree_.empty()) {
        const double rest = std::max(0.0, 1.0 - absSum);
        q.resize(lifeFree_.size());
        w.resize(lifeFree_.size());
        for (std::size_t k = 0; k < lifeFree_.size(); ++k) q[k] = p[intensityPar(lifeFree_[k])];
        hypersphericalToFractions(q.data(), q.size(), w.data());
        double freePs = 0.0;   // p-Ps accompanying the free o-Ps, relative to the free group
        for (std::size_t k = 0; k < lifeFree_.size(); ++k)
            if (isOPs(lifeFree_[k])) freePs += ratio * w[k];
        for (std::size_t k = 0; k < lifeFree_.size(); ++k)
            f.lifetime[lifeFree_[k]] = rest * w[k] / (1.0 + freePs);
        pPsFromFree = rest * freePs / (1.0 + freePs);
    }
    if (withPs) f.lifetime[pPs_] = pPsFromAbs + pPsFromFree;
}

LifetimeModel::Fractions LifetimeModel::fractions(const double* p) const
{
    Fractions f;
    computeFractions(p, f);
    return f;
}

void LifetimeModel::setInitialIntensityParameters(std::vector<double>& p) const
{
    p.resize(nPar(), 0.0);
    std::vector<double> w, q;

    resLayout_.initialFractionParameters(p);

    for (int j : lifeAbs_) p[intensityPar(j)] = lifetimes_[j].intensity;
    if (pPs_ >= 0) p[intensityPar(pPs_)] = lifetimes_[pPs_].intensity;
    w.clear();
    for (int j : lifeFree_) w.push_back(lifetimes_[j].intensity);
    q.resize(w.size());
    fractionsToHyperspherical(w.data(), w.size(), q.data());
    for (std::size_t k = 0; k < lifeFree_.size(); ++k) p[intensityPar(lifeFree_[k])] = q[k];
}

std::vector<std::string> LifetimeModel::parameterNames() const
{
    std::vector<std::string> names(nPar());
    names[kBackground] = "background";
    names[kArea] = "area";
    auto groupName = [](const std::vector<int>& group, int comp, const std::string& absName) {
        for (std::size_t k = 0; k < group.size(); ++k) {
            if (group[k] != comp) continue;
            if (k + 1 == group.size()) return std::string("unused");
            return k == 0 ? std::string("r") : "phi" + std::to_string(k);
        }
        return absName;
    };
    for (int i = 0; i < nResolution(); ++i) {
        const std::string s = std::to_string(i + 1);
        names[sigmaPar(i)] = "sigma_" + s;
        names[offsetPar(i)] = "t0_" + s;
        std::vector<int> resFree;
        for (int k = 0; k < nResolution(); ++k)
            if (resolution_[k].intensityMode == ParMode::Free) resFree.push_back(k);
        names[resolutionFractionPar(i)] = "res_" + s + "_" + groupName(resFree, i, "fraction");
    }
    for (int j = 0; j < nLifetimes(); ++j) {
        const std::string s = std::to_string(j + 1);
        names[tauPar(j)] = "tau_" + s;
        names[intensityPar(j)] = "int_" + s + "_" +
            (lifetimes_[j].pPs ? std::string("pPs_to_oPs_ratio") : groupName(lifeFree_, j, "intensity"));
        if (widthPar_[j] >= 0) names[widthPar_[j]] = "width_" + s;
    }
    return names;
}

double LifetimeModel::operator()(const double* x, const double* p) const
{
    const int n = nPar();
    if (cachedPars_.size() != static_cast<std::size_t>(n) || !std::equal(p, p + n, cachedPars_.begin())) {
        cachedPars_.assign(p, p + n);
        computeFractions(p, cached_);
    }
    const double t = x[0];
    double sum = 0.0;
    for (int i = 0; i < nResolution(); ++i) {
        const double wi = cached_.resolution[i];
        if (wi == 0.0) continue;
        const double sigma = p[sigmaPar(i)];
        const double t0 = p[offsetPar(i)];
        double inner = 0.0;
        for (int j = 0; j < nLifetimes(); ++j) {
            if (componentMask_ >= 0 && j != componentMask_) continue;
            const double Ij = cached_.lifetime[j];
            if (Ij == 0.0) continue;
            const double w = widthPar_[j] >= 0 ? p[widthPar_[j]] : 0.0;
            if (w <= 1e-6) {
                inner += Ij * emg(t, sigma, t0, p[tauPar(j)]);
            } else {
                // log-normal in tau with mean tau_j: ln tau ~ N(ln tau_j - w^2/2, w^2)
                static const GaussHermite gh(32);
                const double mu = std::log(p[tauPar(j)]) - 0.5 * w * w;
                double avg = 0.0;
                for (std::size_t q = 0; q < gh.x.size(); ++q) avg += gh.w[q] * emg(t, sigma, t0, std::exp(mu + w * gh.x[q]));
                inner += Ij * avg;
            }
        }
        sum += wi * inner;
    }
    return (includeBackground_ ? p[kBackground] : 0.0) + p[kArea] * binWidth_ * sum;
}

} // namespace pals
