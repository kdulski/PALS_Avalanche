#include "FitSession.h"
#include "Fitter.h"
#include "FormulaModel.h"
#include "LifetimeModel.h"
#include "ModelSetup.h"
#include "RootFormula.h"
#include "RootIO.h"

#include <TF1.h>
#include <TH1D.h>
#include <TObject.h>

#include <ctime>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>

namespace fs = std::filesystem;

namespace pals {

namespace {

std::string now(const char* format)
{
    std::time_t t = std::time(nullptr);
    char buf[64];
    std::strftime(buf, sizeof(buf), format, std::localtime(&t));
    return buf;
}

} // namespace

Spectrum loadSpectrum(const std::string& path, const Config& cfg)
{
    std::cout << "------------- Reading data -------------\n";
    switch (cfg.input.format) {
        case InputFormat::Root:
            return readRootHistogram(path, cfg.input.rootObject, cfg.histogram.binWidth);
        case InputFormat::Histogram:
            return spectrumFromCounts(readTextColumn(path, cfg.input.column), cfg.histogram);
        case InputFormat::Times:
        default:
            return binTimes(readTextColumn(path, cfg.input.column), cfg.histogram);
    }
}

FitSession::FitSession() = default;
FitSession::~FitSession() = default;

void FitSession::resetFit()
{
    extras_.clear();
    outcome_.reset();
    fitter_.reset();
    lifetimeModel_.reset();
    formulaModel_.reset();
    hasDist_ = false;
}

void FitSession::setConfig(const Config& cfg, const std::string& path)
{
    resetFit();
    est_.reset();
    cfg_ = cfg;
    cfgPath_ = path;
}

void FitSession::loadData(const std::string& path)
{
    resetFit();
    est_.reset();
    auto s = std::make_unique<Spectrum>(loadSpectrum(path, cfg_));
    spectrum_ = std::move(s);
    dataPath_ = path;
}

const Estimates& FitSession::estimate()
{
    if (!spectrum_) throw std::runtime_error("No data loaded");
    resetFit();
    auto e = std::make_unique<Estimates>(pals::estimate(*spectrum_, cfg_.estimation));
    est_ = std::move(e);
    return *est_;
}

const FitOutcome& FitSession::fit(std::function<void(const StageSummary&)> onStage)
{
    estimate();                                   // also discards a previous fit
    const Spectrum& s = *spectrum_;
    const Estimates& est = *est_;
    const auto resolution = resolvedResolution(cfg_, est);

    if (cfg_.model.type == ModelType::Lifetimes) {
        lifetimeModel_ = std::make_unique<LifetimeModel>(resolution, cfg_.lifetimes, s.binWidth, cfg_.fit.oPsThreshold);
        fitter_ = std::make_unique<Fitter>(cfg_, s, est, *lifetimeModel_, lifetimeSpecs(*lifetimeModel_, cfg_, est));
        if (cfg_.fit.profileWidths) {
            std::vector<int> widths;
            for (int j = 0; j < lifetimeModel_->nLifetimes(); ++j)
                if (lifetimeModel_->widthPar(j) >= 0) widths.push_back(lifetimeModel_->widthPar(j));
            fitter_->setProfileParameters(widths);
        }
    } else {
        std::vector<std::string> parNames;
        FormulaFunction formula = compileFormula(cfg_.model.formula, parNames);
        const bool convolve = cfg_.model.convolve;
        const double step = convolve ? formulaGridStep(cfg_, resolution) : 0.0;
        const double range = s.center(est.fitTo) - est.peakPosition;
        formulaModel_ = std::make_unique<FormulaModel>(formula, parNames, convolve,
                                                       convolve ? resolution : std::vector<ResolutionComponent>{},
                                                       s.binWidth, step, cfg_.model.windowSigma, range);
        if (cfg_.distribution.enabled)
            std::cerr << "Warning: the distribution analysis is available only for the lifetime model (ignored)\n";
        fitter_ = std::make_unique<Fitter>(cfg_, s, est, *formulaModel_, formulaSpecs(*formulaModel_, cfg_, s, est));
    }
    fitter_->setStageCallback(std::move(onStage));
    outcome_ = std::make_unique<FitOutcome>(fitter_->run());
    fitDate_ = now("%Y-%m-%d %H:%M:%S");

    if (lifetimeModel_) {
        if (cfg_.distribution.enabled) {
            std::cout << "\n------------- Lifetime distribution -------------\n";
            dist_ = analyzeDistribution(distributionProblem(*lifetimeModel_, *outcome_, s, est), cfg_.distribution,
                                        distributionRegions(cfg_));
            hasDist_ = true;
        }
        extras_ = lifetimeComponentFunctions(*lifetimeModel_, *outcome_, fitter_->function()->GetNpx());
        if (hasDist_)
            for (auto& h : distributionHistograms(dist_)) extras_.push_back(std::move(h));
    }
    return *outcome_;
}

TH1D* FitSession::histogram() const { return fitter_ ? fitter_->histogram() : nullptr; }
TF1* FitSession::function() const { return fitter_ ? fitter_->function() : nullptr; }

std::vector<TF1*> FitSession::components() const
{
    std::vector<TF1*> out;
    for (const auto& o : extras_)
        if (auto* f = dynamic_cast<TF1*>(o.get())) out.push_back(f);
    return out;
}

void FitSession::writeReport(std::ostream& os) const
{
    if (!outcome_) throw std::runtime_error("No fit to report");
    const ReportInfo info{dataPath_, cfgPath_.empty() ? std::string("(from the GUI)") : cfgPath_, fitDate_};
    const FitOutcome& fit = *outcome_;
    const ResultsSection results = [&](std::ostream& o) {
        if (lifetimeModel_) {
            writeLifetimeResults(o, *lifetimeModel_, *est_, fit);
            if (hasDist_) writeDistributionReport(o, dist_, cfg_.distribution, fit.params[LifetimeModel::kArea]);
        } else {
            writeFormulaResults(o, *formulaModel_, cfg_.model.formula, *est_, fit);
        }
    };
    pals::writeReport(os, info, *spectrum_, *est_, fit, results);
}

void FitSession::save(const std::string& basePath) const
{
    if (!outcome_) throw std::runtime_error("No fit to save");
    const fs::path parent = fs::path(basePath).parent_path();
    if (!parent.empty()) fs::create_directories(parent);
    std::ofstream txt(basePath + ".txt");
    if (!txt) throw std::runtime_error("Cannot write " + basePath + ".txt");
    writeReport(txt);
    // the configuration actually used (also when it was edited in the GUI)
    writeRootFile(basePath + ".root", *fitter_, *outcome_, cfg_.toJson(), extras_);
}

std::string FitSession::defaultOutputBase() const
{
    std::string stem = dataPath_.empty() ? std::string("spectrum") : fs::path(dataPath_).stem().string();
    if (cfg_.output.timestamp) stem += "_" + now("%Y%m%d-%H%M%S");
    return (fs::path(cfg_.output.directory) / stem).string();
}

} // namespace pals
