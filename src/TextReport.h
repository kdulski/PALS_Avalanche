#ifndef PALS_TEXTREPORT_H
#define PALS_TEXTREPORT_H

#include "FormulaModel.h"
#include "LifetimeModel.h"
#include "Spectrum.h"

#include <functional>
#include <iosfwd>
#include <string>
#include <vector>

namespace pals {

struct StageSummary {
    int stage = 0;
    std::string options;
    int status = 0;
    double chi2 = 0.0;
    int ndf = 0;
};

// 68 % profile-likelihood interval of one parameter (objective rises by 1).
struct ProfileInterval {
    int par = -1;
    double best = 0.0, lower = 0.0, upper = 0.0;
    bool lowerAtLimit = false, upperAtLimit = false;
    bool lowerMinimumFound = false;      // a better minimum was met while profiling
    int fits = 0;
};

// Everything needed to report a fit, independent of ROOT.
struct FitOutcome {
    std::vector<std::string> names;
    std::vector<double> params, errors, lower, upper;
    std::vector<bool> fixed;
    std::vector<std::vector<double>> covariance;   // nPar x nPar, zero rows for fixed parameters
    double chi2 = 0.0;
    int ndf = 0;
    double pValue = 0.0;
    int status = -1;
    int covarianceStatus = -1;
    double fitMin = 0.0, fitMax = 0.0;
    int fitBins = 0;
    std::vector<StageSummary> stages;
    std::vector<ProfileInterval> profiles;
    std::string objective;               // statistic used for the profiles
};

struct ValueError {
    double value = 0.0;
    double error = 0.0;
};

// Values g(p) with errors from the full covariance matrix: sigma^2 = J C J^T
// (numerical Jacobian over the free parameters).
std::vector<ValueError> propagate(const std::function<std::vector<double>(const double*)>& g, const FitOutcome& fit);

// Physical fractions of the lifetime model (resolution weights, lifetime intensities).
struct DerivedFractions {
    std::vector<ValueError> resolution;
    std::vector<ValueError> lifetime;
};
DerivedFractions derivedFractions(const LifetimeModel& model, const FitOutcome& fit);

struct ReportInfo {
    std::string dataPath;
    std::string configPath;
    std::string date;
};

// Model-specific "Results" section
using ResultsSection = std::function<void(std::ostream&)>;

void writeReport(std::ostream& os, const ReportInfo& info, const Spectrum& spectrum, const Estimates& est,
                 const FitOutcome& fit, const ResultsSection& results);

// Note on the statistics needed by width / distribution analyses (from simulations).
std::string statisticsNote(double counts, double low, double good, const std::string& what);

void writeLifetimeResults(std::ostream& os, const LifetimeModel& model, const Estimates& est, const FitOutcome& fit);
void writeFormulaResults(std::ostream& os, const FormulaModel& model, const std::string& formula,
                         const Estimates& est, const FitOutcome& fit);

} // namespace pals

#endif
