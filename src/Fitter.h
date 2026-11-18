#ifndef PALS_FITTER_H
#define PALS_FITTER_H

#include "Config.h"
#include "FitParameters.h"
#include "Spectrum.h"
#include "TextReport.h"

#include <TF1.h>
#include <TFitResultPtr.h>

#include <functional>
#include <memory>
#include <vector>

class TH1D;

namespace pals {

// Staged fit of any model (functor usable by TF1) described by ParameterSpecs.
//
//   stage 0           : fit.first_pass_options
//   stage 1..N        : fit.options, N = fit.iterations
// In every stage the parameters are fixed / released according to their
// ParameterSpec (mode, limits, release stage); "partial" parameters get a
// window of +- partial_variation around their current value.
class Fitter {
public:
    template <class Model>
    Fitter(const Config& cfg, const Spectrum& spectrum, const Estimates& est, const Model& model,
           std::vector<ParameterSpec> specs)
        : Fitter(cfg, spectrum, est, std::move(specs))
    {
        func_ = std::make_unique<TF1>("fit", model, fitMin(), fitMax(), static_cast<int>(specs_.size()), 1,
                                      TF1::EAddToList::kNo);
        initialiseFunction();
    }
    ~Fitter();

    FitOutcome run();

    // parameters for which run() computes profile-likelihood intervals
    void setProfileParameters(std::vector<int> pars) { profilePars_ = std::move(pars); }
    // called after every fit stage (progress display)
    void setStageCallback(std::function<void(const StageSummary&)> cb) { onStage_ = std::move(cb); }

    TH1D* histogram() const { return hist_.get(); }
    TF1* function() const { return func_.get(); }
    double fitMin() const;
    double fitMax() const;

private:
    Fitter(const Config& cfg, const Spectrum& spectrum, const Estimates& est, std::vector<ParameterSpec> specs);
    void initialiseFunction();
    void configureStage(int stage);
    void fix(int ipar);
    void limit(int ipar, double lo, double hi);
    void partial(int ipar);
    StageSummary fitStage(int stage, const std::string& options);
    FitOutcome collect() const;
    double objective() const;
    std::string objectiveName() const;
    void applyState(const FitOutcome& state);
    ProfileInterval profile(int ipar, const FitOutcome& base);

    const Config& cfg_;
    const Spectrum& spectrum_;
    const Estimates& est_;
    std::vector<ParameterSpec> specs_;
    std::unique_ptr<TH1D> hist_;
    std::unique_ptr<TF1> func_;
    TFitResultPtr result_;
    std::vector<StageSummary> stages_;
    std::vector<int> profilePars_;
    std::function<void(const StageSummary&)> onStage_;
};

} // namespace pals

#endif
