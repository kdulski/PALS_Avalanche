#ifndef PALS_MODELSETUP_H
#define PALS_MODELSETUP_H

// Initial values, limits and release stages of the fit parameters for the
// lifetime model and for the formula model (ROOT-independent).

#include "Config.h"
#include "Distribution.h"
#include "FitParameters.h"
#include "FormulaModel.h"
#include "LifetimeModel.h"
#include "Spectrum.h"
#include "TextReport.h"

#include <string>
#include <vector>

namespace pals {

// Resolution components with sigma = "auto" replaced by the estimate.
std::vector<ResolutionComponent> resolvedResolution(const Config& cfg, const Estimates& est);

// Values available as "estimate:<name>" in formula parameters:
//   background, background_error, area, peak_position, peak_counts, sigma,
//   bin_width, fit_start, fit_end
double estimateValue(const std::string& name, const Spectrum& s, const Estimates& est);

std::vector<ParameterSpec> lifetimeSpecs(const LifetimeModel& model, const Config& cfg, const Estimates& est);

// Grid step of the numerical convolution: model.step or min(sigma)/points_per_sigma.
double formulaGridStep(const Config& cfg, const std::vector<ResolutionComponent>& resolution);

std::vector<ParameterSpec> formulaSpecs(const FormulaModel& model, const Config& cfg, const Spectrum& s,
                                        const Estimates& est);

// Grid problem for the distribution analysis from the result of the lifetime-model fit
// (same fit range, resolution, t0, background; Pearson weights from the fitted model).
DistributionProblem distributionProblem(const LifetimeModel& model, const FitOutcome& fit, const Spectrum& s,
                                        const Estimates& est);

// Regions of the distribution summaries: from the config or [o-Ps threshold, grid max].
std::vector<std::pair<double, double>> distributionRegions(const Config& cfg);

} // namespace pals

#endif
