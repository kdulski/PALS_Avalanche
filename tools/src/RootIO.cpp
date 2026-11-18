#include "RootIO.h"

#include <TF1.h>
#include <TFile.h>
#include <TH1.h>
#include <TH1D.h>
#include <TNamed.h>
#include <TObject.h>

#include <cmath>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>

namespace pals {

Spectrum readRootHistogram(const std::string& path, const std::string& object, double binWidth)
{
    std::unique_ptr<TFile> file(TFile::Open(path.c_str(), "READ"));
    if (!file || file->IsZombie()) throw std::runtime_error("Cannot open ROOT file " + path);
    auto* h = dynamic_cast<TH1*>(file->Get(object.c_str()));
    if (!h) throw std::runtime_error("No TH1 \"" + object + "\" in " + path);
    if (h->GetDimension() != 1) throw std::runtime_error("\"" + object + "\" is not a 1D histogram");
    if (h->GetXaxis()->IsVariableBinSize()) throw std::runtime_error("Variable bin widths are not supported");

    std::unique_ptr<TH1> hist(static_cast<TH1*>(h->Clone("input")));
    hist->SetDirectory(nullptr);

    const double original = hist->GetXaxis()->GetBinWidth(1);
    const double ratio = binWidth / original;
    const long rebin = std::lround(ratio);
    if (rebin < 1 || std::fabs(ratio - rebin) > 1e-3 * ratio) {
        std::ostringstream msg;
        msg << "bin_width " << binWidth << " ns is not an integer multiple of the histogram bin width "
            << original << " ns";
        throw std::runtime_error(msg.str());
    }
    if (rebin > 1) {
        std::cout << "Rebinning " << object << " by " << rebin << "\n";
        hist->Rebin(static_cast<int>(rebin));
    }

    Spectrum s;
    s.min = hist->GetXaxis()->GetXmin();
    s.binWidth = hist->GetXaxis()->GetBinWidth(1);
    const int n = hist->GetNbinsX();
    s.counts.resize(n);
    for (int i = 1; i <= n; ++i) s.counts[i - 1] = hist->GetBinContent(i);   // under/overflow skipped
    std::cout << "Histogram " << object << ": " << n << " bins of " << s.binWidth << " ns in [" << s.min
              << ", " << s.max() << "] ns\n";
    return s;
}

std::vector<std::unique_ptr<TObject>> lifetimeComponentFunctions(const LifetimeModel& model, const FitOutcome& fit, int npx)
{
    const int colors[] = {kBlue, kGreen + 2, kMagenta, kOrange + 7, kCyan + 2, kViolet};
    std::vector<std::unique_ptr<TObject>> out;
    for (int j = 0; j < model.nLifetimes(); ++j) {
        LifetimeModel single = model;
        single.setComponentMask(j);
        single.setIncludeBackground(false);
        const std::string name = "component_" + std::to_string(j + 1);
        auto f = std::make_unique<TF1>(name.c_str(), single, fit.fitMin, fit.fitMax, model.nPar(), 1,
                                       TF1::EAddToList::kNo);
        f->SetParameters(fit.params.data());
        f->SetNpx(npx);
        f->SetLineColor(colors[j % 6]);
        f->SetLineStyle(2);
        out.push_back(std::move(f));
    }
    return out;
}

std::vector<std::unique_ptr<TObject>> distributionHistograms(const DistributionAnalysis& a)
{
    const auto& g = a.grid;
    const std::size_t n = g.size();
    std::vector<double> edges(n + 1);
    const double r = std::sqrt(g[1] / g[0]);
    for (std::size_t k = 0; k < n; ++k) edges[k] = g[k] / r;
    edges[n] = g[n - 1] * r;
    std::vector<std::unique_ptr<TObject>> out;
    auto hist = [&](const std::string& name, const std::string& title, const std::vector<double>& d) {
        auto h = std::make_unique<TH1D>(name.c_str(), (title + ";lifetime [ns];density per ln#tau").c_str(),
                                        static_cast<int>(n), edges.data());
        h->SetDirectory(nullptr);
        for (std::size_t k = 0; k < n; ++k) h->SetBinContent(static_cast<int>(k) + 1, d[k]);
        out.push_back(std::move(h));
    };
    for (const auto& e : a.estimates) {
        const std::string key = e.method == "MaxEnt" ? "maxent" : "nnls_bagged";
        hist("distribution_" + key, e.method, densityPerLnTau(g, e.x));
        if (!e.densityLow.empty()) {
            hist("distribution_" + key + "_q16", e.method + " bootstrap 16 %", e.densityLow);
            hist("distribution_" + key + "_q84", e.method + " bootstrap 84 %", e.densityHigh);
        }
    }
    return out;
}

void writeRootFile(const std::string& path, const Fitter& fitter, const FitOutcome& fit,
                   const std::string& configText, const std::vector<std::unique_ptr<TObject>>& extra)
{
    std::unique_ptr<TFile> file(TFile::Open(path.c_str(), "RECREATE"));
    if (!file || file->IsZombie()) throw std::runtime_error("Cannot create ROOT file " + path);

    TH1D* data = fitter.histogram();
    TF1* total = fitter.function();
    data->Write("data");
    total->SetLineColor(kRed);
    total->Write("fit");
    for (const auto& f : extra) f->Write();

    TF1 bg("background", "[0]", fit.fitMin, fit.fitMax, TF1::EAddToList::kNo);
    bg.SetParameter(0, fit.params[0]);
    bg.SetLineColor(kGray + 2);
    bg.SetLineStyle(3);
    bg.Write();

    // normalised residuals (data - fit)/sqrt(data) in the fit range
    TH1D residuals("residuals", "Residuals;time [ns];(data - fit) / #sqrt{data}", data->GetNbinsX(),
                   data->GetXaxis()->GetXmin(), data->GetXaxis()->GetXmax());
    residuals.SetDirectory(nullptr);
    for (int i = 1; i <= data->GetNbinsX(); ++i) {
        const double x = data->GetBinCenter(i);
        const double y = data->GetBinContent(i);
        if (x < fit.fitMin || x > fit.fitMax || y <= 0.0) continue;
        residuals.SetBinContent(i, (y - total->Eval(x)) / std::sqrt(y));
    }
    residuals.Write();

    TNamed("config", configText.c_str()).Write();
    file->Close();
}

} // namespace pals
