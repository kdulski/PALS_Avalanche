#ifndef PALS_RESOLUTION_H
#define PALS_RESOLUTION_H

#include "Components.h"
#include "Intensities.h"

#include <vector>

namespace pals {

// Resolution function = sum of Gaussians; where their parameters live in the
// parameter vector of a model.
struct ResolutionLayout {
    std::vector<ResolutionComponent> components;
    std::vector<int> sigmaPar, offsetPar, fractionPar;

    int size() const { return static_cast<int>(components.size()); }

    // fraction parameters of the hyperspherical group (free intensities), in order
    std::vector<int> groupPars() const
    {
        std::vector<int> g;
        for (int i = 0; i < size(); ++i)
            if (components[i].intensityMode == ParMode::Free) g.push_back(fractionPar[i]);
        return g;
    }

    std::vector<double> weights(const double* p) const
    {
        std::vector<double> pars(size());
        std::vector<bool> isFree(size());
        for (int i = 0; i < size(); ++i) {
            pars[i] = p[fractionPar[i]];
            isFree[i] = components[i].intensityMode == ParMode::Free;
        }
        return groupFractions(pars, isFree);
    }

    // initial values of the fraction parameters from the component intensities
    void initialFractionParameters(std::vector<double>& p) const
    {
        std::vector<double> w;
        for (int i = 0; i < size(); ++i)
            if (components[i].intensityMode == ParMode::Free) w.push_back(components[i].intensity);
            else p[fractionPar[i]] = components[i].intensity;
        std::vector<double> q(w.size());
        fractionsToHyperspherical(w.data(), w.size(), q.data());
        const auto g = groupPars();
        for (std::size_t k = 0; k < g.size(); ++k) p[g[k]] = q[k];
    }
};

} // namespace pals

#endif
