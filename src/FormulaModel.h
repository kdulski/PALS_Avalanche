#ifndef PALS_FORMULAMODEL_H
#define PALS_FORMULAMODEL_H

#include "Convolution.h"
#include "Resolution.h"

#include <functional>
#include <string>
#include <vector>

namespace pals {

// D(x; params): the user formula (in the program a compiled TFormula).
using FormulaFunction = std::function<double(double x, const double* params)>;

// Model given by a formula from the config.
//
//   convolve = true:   N(t) = B + binWidth * sum_i w_i  int_0^inf D(x) G(t - t0_i - x; sigma_i) dx
//                      x = time since t0_i, D = 0 for x < 0
//   convolve = false:  N(t) = B + binWidth * D(t)
//                      (resolution is part of the formula, e.g. via pals_emg)
//
// D is a density in counts per ns (multiplied by the bin width it gives counts
// per bin), so e.g. D = A*I/tau*exp(-x/tau) has A*I counts.
//
// Parameter layout:
//   0                  background B
//   1 + 3i, 2 + 3i,    sigma_i, t0_i, fraction parameter of Gaussian i
//   3 + 3i             (only with convolve = true)
//   first + k          k-th parameter of the formula
class FormulaModel {
public:
    FormulaModel(FormulaFunction formula, std::vector<std::string> formulaParNames, bool convolve,
                 std::vector<ResolutionComponent> resolution, double binWidth,
                 double step, double windowSigma, double initialRange);

    static constexpr int kBackground = 0;

    bool convolve() const { return convolve_; }
    int nResolution() const { return resLayout_.size(); }
    int nFormulaPars() const { return static_cast<int>(names_.size()); }
    int firstFormulaPar() const { return 1 + 3 * nResolution(); }
    int nPar() const { return firstFormulaPar() + nFormulaPars(); }
    int sigmaPar(int i) const { return 1 + 3 * i; }
    int offsetPar(int i) const { return 2 + 3 * i; }
    int resolutionFractionPar(int i) const { return 3 + 3 * i; }
    const ResolutionLayout& resolutionLayout() const { return resLayout_; }
    const std::vector<std::string>& formulaParNames() const { return names_; }
    double step() const { return step_; }

    std::vector<std::string> parameterNames() const;

    double operator()(const double* x, const double* p) const;

private:
    void prepare(const double* formulaPars, double neededRange) const;

    FormulaFunction formula_;
    std::vector<std::string> names_;
    bool convolve_;
    ResolutionLayout resLayout_;
    double binWidth_;
    double step_;
    double window_;
    double initialRange_;

    // D on the grid, recomputed only when the formula parameters change
    mutable GaussianConvolution grid_;
    mutable std::vector<double> cachedPars_;
};

} // namespace pals

#endif
