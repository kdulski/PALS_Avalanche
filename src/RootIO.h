#ifndef PALS_ROOTIO_H
#define PALS_ROOTIO_H

#include "Fitter.h"
#include "Distribution.h"
#include "LifetimeModel.h"

#include <TObject.h>
#include "Spectrum.h"
#include "TextReport.h"

#include <memory>
#include <string>
#include <vector>

namespace pals {

// Reads a TH1 from a ROOT file (object path like "dir/histogram") and rebins
// it to binWidth (must be an integer multiple of the original bin width).
Spectrum readRootHistogram(const std::string& path, const std::string& object, double binWidth);



// Single lifetime components (resolution-convolved, without background) for drawing.
std::vector<std::unique_ptr<TObject>> lifetimeComponentFunctions(const LifetimeModel& model, const FitOutcome& fit, int npx);

// Distributions (density per ln tau, unit total intensity) as histograms with
// log-spaced bins in tau; bagged NNLS also with its 16 / 84 % bootstrap band.
std::vector<std::unique_ptr<TObject>> distributionHistograms(const DistributionAnalysis& a);

// Writes the data histogram, the fitted function, the background, the
// normalised residuals, the config text and the extra functions (components).
void writeRootFile(const std::string& path, const Fitter& fitter, const FitOutcome& fit,
                   const std::string& configText, const std::vector<std::unique_ptr<TObject>>& extra = {});

} // namespace pals

#endif
