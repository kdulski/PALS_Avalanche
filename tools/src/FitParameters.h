#ifndef PALS_FITPARAMETERS_H
#define PALS_FITPARAMETERS_H

#include "Components.h"

#include <string>
#include <vector>

namespace pals {

// Description of one fit parameter for the staged fit (independent of the model).
//
// In stage s the parameter is
//   fixed            if mode == Fixed, or s < releaseStage,
//                    or mode == Partial and s == 0
//   in a window      if mode == Partial: value +- partial_variation, re-centred each stage
//   within limits    otherwise (or unbounded when hasLimits == false)
struct ParameterSpec {
    std::string name;
    double value = 0.0;
    bool hasLimits = false;
    double lower = 0.0, upper = 0.0;
    ParMode mode = ParMode::Free;
    int releaseStage = 0;
};

inline ParameterSpec makeSpec(std::string name, double value, ParMode mode = ParMode::Free, int releaseStage = 0)
{
    ParameterSpec s;
    s.name = std::move(name);
    s.value = value;
    s.mode = mode;
    s.releaseStage = releaseStage;
    return s;
}

inline ParameterSpec makeSpec(std::string name, double value, double lower, double upper,
                              ParMode mode = ParMode::Free, int releaseStage = 0)
{
    ParameterSpec s = makeSpec(std::move(name), value, mode, releaseStage);
    s.hasLimits = true;
    s.lower = lower;
    s.upper = upper;
    return s;
}

} // namespace pals

#endif
