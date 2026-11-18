#ifndef PALS_FITSESSION_H
#define PALS_FITSESSION_H

// One analysis: configuration, data, preliminary estimates, fit, report and output
// files. Used by the command-line program and by the GUI.

#include "Config.h"
#include "Distribution.h"
#include "Spectrum.h"
#include "TextReport.h"

#include <functional>
#include <iosfwd>
#include <memory>
#include <string>
#include <vector>

class TF1;
class TH1D;
class TObject;

namespace pals {

class Fitter;
class FormulaModel;
class LifetimeModel;

// Reads a spectrum according to the input / histogram part of the configuration.
Spectrum loadSpectrum(const std::string& path, const Config& cfg);

class FitSession {
public:
    FitSession();
    ~FitSession();

    // A new configuration discards the estimates and the fit (not the data).
    void setConfig(const Config& cfg, const std::string& path = "");
    const Config& config() const { return cfg_; }
    const std::string& configPath() const { return cfgPath_; }

    void loadData(const std::string& path);      // with the current input / histogram settings
    bool hasData() const { return static_cast<bool>(spectrum_); }
    const std::string& dataPath() const { return dataPath_; }
    const Spectrum& spectrum() const { return *spectrum_; }

    const Estimates& estimate();                 // from the current config; discards the fit
    bool hasEstimates() const { return static_cast<bool>(est_); }
    const Estimates& estimates() const { return *est_; }

    const FitOutcome& fit(std::function<void(const StageSummary&)> onStage = {});
    bool hasFit() const { return static_cast<bool>(outcome_); }
    const FitOutcome& outcome() const { return *outcome_; }

    TH1D* histogram() const;                     // data with the fit attached (after a fit)
    TF1* function() const;                       // fitted function (after a fit)
    std::vector<TF1*> components() const;        // single lifetime components (after a fit)

    void writeReport(std::ostream& os) const;
    void save(const std::string& basePath) const;   // basePath.txt and basePath.root
    std::string defaultOutputBase() const;          // <output dir>/<data stem>[_<timestamp>]

private:
    void resetFit();

    Config cfg_;
    std::string cfgPath_, dataPath_;
    std::unique_ptr<Spectrum> spectrum_;
    std::unique_ptr<Estimates> est_;
    std::unique_ptr<LifetimeModel> lifetimeModel_;
    std::unique_ptr<FormulaModel> formulaModel_;
    std::unique_ptr<Fitter> fitter_;
    std::unique_ptr<FitOutcome> outcome_;
    DistributionAnalysis dist_;
    bool hasDist_ = false;
    std::vector<std::unique_ptr<TObject>> extras_;
    std::string fitDate_;
};

} // namespace pals

#endif
