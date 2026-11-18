#include "FormulaModel.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace pals {

FormulaModel::FormulaModel(FormulaFunction formula, std::vector<std::string> formulaParNames, bool convolve,
                           std::vector<ResolutionComponent> resolution, double binWidth,
                           double step, double windowSigma, double initialRange)
    : formula_(std::move(formula)), names_(std::move(formulaParNames)), convolve_(convolve),
      binWidth_(binWidth), step_(step), window_(windowSigma), initialRange_(initialRange)
{
    if (convolve_) {
        if (resolution.empty()) throw std::invalid_argument("FormulaModel: convolution needs resolution components");
        if (step_ <= 0.0) throw std::invalid_argument("FormulaModel: grid step must be > 0");
        resLayout_.components = std::move(resolution);
        for (int i = 0; i < static_cast<int>(resLayout_.components.size()); ++i) {
            resLayout_.sigmaPar.push_back(sigmaPar(i));
            resLayout_.offsetPar.push_back(offsetPar(i));
            resLayout_.fractionPar.push_back(resolutionFractionPar(i));
        }
    }
}

std::vector<std::string> FormulaModel::parameterNames() const
{
    std::vector<std::string> n(nPar());
    n[kBackground] = "background";
    const auto group = resLayout_.groupPars();
    for (int i = 0; i < nResolution(); ++i) {
        const std::string s = std::to_string(i + 1);
        n[sigmaPar(i)] = "sigma_" + s;
        n[offsetPar(i)] = "t0_" + s;
        std::string role = "fraction";
        for (std::size_t k = 0; k < group.size(); ++k)
            if (group[k] == resolutionFractionPar(i))
                role = k + 1 == group.size() ? "unused" : (k == 0 ? "r" : "phi" + std::to_string(k));
        n[resolutionFractionPar(i)] = "res_" + s + "_" + role;
    }
    for (int k = 0; k < nFormulaPars(); ++k) n[firstFormulaPar() + k] = names_[k];
    return n;
}

void FormulaModel::prepare(const double* fp, double neededRange) const
{
    const auto n = static_cast<std::size_t>(nFormulaPars());
    const bool same = cachedPars_.size() == n && std::equal(fp, fp + n, cachedPars_.begin()) && grid_.size() > 0;
    if (!same) {
        cachedPars_.assign(fp, fp + n);
        grid_.reset(step_);
    }
    const double target = std::max(neededRange, same ? 0.0 : initialRange_);
    if (grid_.size() > 0 && grid_.range() >= target) return;
    // extend the grid (odd number of nodes, i.e. whole panels); 25% margin
    const double h = step_;
    const double goal = std::max(target, 1.25 * grid_.range());
    if (grid_.size() == 0) grid_.append(formula_(0.0, fp));
    while (grid_.range() < goal) {
        const std::size_t k = grid_.size();
        grid_.append(formula_(k * h, fp));
        grid_.append(formula_((k + 1) * h, fp));
    }
}

double FormulaModel::operator()(const double* x, const double* p) const
{
    const double t = x[0];
    const double* fp = p + firstFormulaPar();
    if (!convolve_) return p[kBackground] + binWidth_ * formula_(t, fp);

    double needed = 0.0;
    for (int i = 0; i < nResolution(); ++i)
        needed = std::max(needed, t - p[offsetPar(i)] + window_ * p[sigmaPar(i)]);
    prepare(fp, needed);

    const auto w = resLayout_.weights(p);
    double sum = 0.0;
    for (int i = 0; i < nResolution(); ++i) {
        if (w[i] == 0.0) continue;
        sum += w[i] * grid_(t - p[offsetPar(i)], p[sigmaPar(i)], window_);
    }
    return p[kBackground] + binWidth_ * sum;
}

} // namespace pals
