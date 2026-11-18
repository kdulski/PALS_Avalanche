#include "Config.h"

#include <nlohmann/json.hpp>

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>

namespace pals {

using json = nlohmann::json;

namespace {

[[noreturn]] void fail(const std::string& msg)
{
    throw std::runtime_error("Config: " + msg);
}

// Warn about keys that are not recognised (typos would otherwise be silently ignored).
void checkKeys(const json& j, const std::string& where, const std::set<std::string>& allowed)
{
    for (auto it = j.begin(); it != j.end(); ++it)
        if (!allowed.count(it.key()))
            std::cerr << "Config warning: unknown key \"" << it.key() << "\" in " << where << " (ignored)\n";
}

const json& section(const json& j, const char* name)
{
    static const json empty = json::object();
    if (!j.contains(name)) return empty;
    if (!j.at(name).is_object()) fail(std::string("\"") + name + "\" must be an object");
    return j.at(name);
}

ParMode parseMode(const std::string& s, const std::string& where)
{
    if (s == "free") return ParMode::Free;
    if (s == "fixed") return ParMode::Fixed;
    if (s == "partial") return ParMode::Partial;
    fail("unknown mode \"" + s + "\" in " + where + " (use free / fixed / partial)");
}

// Legacy component types of FitDetails:
//   nf  - all free             f   - lifetime and intensity fixed
//   pf  - lifetime partially free (intensity free)
//   lpf - lifetime and intensity partially free
//   lf  - lifetime fixed (intensity free)
//   lff - lifetime fixed, intensity partially free
//   ps  - p-Ps, intensity = fraction of the o-Ps intensities
void applyLegacyType(const std::string& type, LifetimeComponent& c, const std::string& where)
{
    using M = ParMode;
    if (type == "nf")       { c.tauMode = M::Free;    c.intensityMode = M::Free; }
    else if (type == "f")   { c.tauMode = M::Fixed;   c.intensityMode = M::Fixed; }
    else if (type == "pf")  { c.tauMode = M::Partial; c.intensityMode = M::Free; }
    else if (type == "lpf") { c.tauMode = M::Partial; c.intensityMode = M::Partial; }
    else if (type == "lf")  { c.tauMode = M::Fixed;   c.intensityMode = M::Free; }
    else if (type == "lff") { c.tauMode = M::Fixed;   c.intensityMode = M::Partial; }
    else if (type == "ps")  { c.tauMode = M::Partial; c.intensityMode = M::Fixed; c.pPs = true; }
    else fail("unknown component type \"" + type + "\" in " + where);
}

template <typename T>
std::pair<T, T> readPair(const json& j, const char* key, std::pair<T, T> def)
{
    if (!j.contains(key)) return def;
    const json& a = j.at(key);
    if (!a.is_array() || a.size() != 2) fail(std::string("\"") + key + "\" must be [low, high]");
    return {a[0].get<T>(), a[1].get<T>()};
}

// Intensities of one group of components (resolution or lifetimes).
// Absolute (fixed/partial) intensities are used as given; the free ones only
// matter relative to each other. Without free components the absolute ones are
// renormalised to 1.
template <typename Comp>
void checkIntensities(std::vector<Comp>& comps, const std::string& what, bool (*skip)(const Comp&))
{
    double absSum = 0.0, freeSum = 0.0;
    int nFree = 0;
    for (const auto& c : comps) {
        if (skip(c)) continue;
        if (c.intensity < 0.0) fail("negative intensity in " + what);
        if (c.intensityMode == ParMode::Free) { freeSum += c.intensity; ++nFree; }
        else absSum += c.intensity;
    }
    if (nFree == 0) {
        if (absSum <= 0.0) fail("all " + what + " intensities are zero");
        if (std::fabs(absSum - 1.0) > 1e-6) {
            std::cerr << "Config warning: " << what << " intensities sum to " << absSum
                      << " and none is free -> renormalised to 1\n";
            for (auto& c : comps) if (!skip(c)) c.intensity /= absSum;
        }
        return;
    }
    if (absSum >= 1.0) fail("fixed/partial " + what + " intensities sum to >= 1, nothing left for the free ones");
    if (freeSum <= 0.0) {
        std::cerr << "Config warning: free " << what << " intensities are zero -> equal shares assumed\n";
        for (auto& c : comps) if (!skip(c) && c.intensityMode == ParMode::Free) c.intensity = 1.0;
    } else if (std::fabs(absSum + freeSum - 1.0) > 0.01) {
        std::cerr << "Config info: " << what << " intensities sum to " << absSum + freeSum
                  << "; free ones are used only as relative shares of " << 1.0 - absSum << "\n";
    }
}

// "model" section of a formula model
void parseFormulaModel(const json& mo, ModelConfig& m)
{
    if (!mo.contains("formula") || !mo.at("formula").is_string() || mo.at("formula").get<std::string>().empty())
        fail("model.formula (string) is required for model.type = formula");
    m.formula = mo.at("formula").get<std::string>();
    m.convolve = mo.value("convolve", true);
    m.pointsPerSigma = mo.value("points_per_sigma", m.pointsPerSigma);
    m.step = mo.value("step", m.step);
    m.windowSigma = mo.value("window_sigma", m.windowSigma);
    if (m.pointsPerSigma <= 0.0 || m.step < 0.0 || m.windowSigma <= 0.0)
        fail("model.points_per_sigma, model.step and model.window_sigma must be positive");

    if (!mo.contains("parameters") || !mo.at("parameters").is_object())
        fail("model.parameters must be an object {\"name\": {...}, ...}");
    for (auto it = mo.at("parameters").begin(); it != mo.at("parameters").end(); ++it) {
        FormulaParameter par;
        par.name = it.key();
        const std::string where = "model.parameters." + par.name;
        const json& d = it.value();
        json value;
        if (d.is_object()) {
            checkKeys(d, where, {"value", "limits", "mode", "stage", "comment"});
            if (!d.contains("value")) fail(where + ".value is required");
            value = d.at("value");
            if (d.contains("limits")) {
                std::tie(par.lower, par.upper) = readPair<double>(d, "limits", {0.0, 0.0});
                if (par.upper <= par.lower) fail(where + ".limits must be [low, high] with low < high");
                par.hasLimits = true;
            }
            if (d.contains("mode")) par.mode = parseMode(d.at("mode").get<std::string>(), where);
            par.stage = d.value("stage", 0);
            if (par.stage < 0) fail(where + ".stage must be >= 0");
        } else {
            value = d;     // shorthand: "name": 0.2
        }
        if (value.is_number()) {
            par.value = value.get<double>();
        } else if (value.is_string()) {
            // "estimate:<name>" or "estimate:<name>*<factor>"
            const std::string v = value.get<std::string>();
            const std::string prefix = "estimate:";
            if (v.compare(0, prefix.size(), prefix) != 0) fail(where + ".value must be a number or \"estimate:<name>\"");
            std::string rest = v.substr(prefix.size());
            const auto star = rest.find('*');
            if (star != std::string::npos) {
                try { par.estimateFactor = std::stod(rest.substr(star + 1)); }
                catch (...) { fail(where + ": cannot read the factor in \"" + v + "\""); }
                rest = rest.substr(0, star);
            }
            par.estimate = rest;
        } else {
            fail(where + ".value must be a number or a string");
        }
        m.parameters.push_back(par);
    }
}

} // namespace

Config Config::load(const std::string& path)
{
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Cannot open config file " + path);
    std::stringstream buffer;
    buffer << in.rdbuf();
    return parse(buffer.str(), path);
}

Config Config::parse(const std::string& text, const std::string& sourcePath)
{
    json root;
    try {
        root = json::parse(text, nullptr, /*allow_exceptions=*/true, /*ignore_comments=*/true);
    } catch (const json::parse_error& e) {
        fail(std::string("JSON syntax error: ") + e.what());
    }
    if (!root.is_object()) fail("top level must be an object");

    Config cfg;
    cfg.sourcePath = sourcePath;
    cfg.sourceText = text;
    if (root.contains("description") && root.at("description").is_string())
        cfg.description = root.at("description").get<std::string>();

    try {
        checkKeys(root, "top level", {"description", "input", "histogram", "estimation", "model",
                                      "resolution", "lifetimes", "fit", "distribution", "output"});

        // --- model type (default: discrete lifetime components)
        const json& mo = section(root, "model");
        checkKeys(mo, "model", {"type", "formula", "convolve", "parameters", "points_per_sigma", "step",
                                "window_sigma"});
        const std::string modelType = mo.value("type", "lifetimes");
        if (modelType == "lifetimes") cfg.model.type = ModelType::Lifetimes;
        else if (modelType == "formula") cfg.model.type = ModelType::Formula;
        else fail("model.type must be lifetimes / formula");
        const bool isFormula = cfg.model.type == ModelType::Formula;
        if (isFormula) parseFormulaModel(mo, cfg.model);
        const bool needResolution = !isFormula || cfg.model.convolve;

        // --- input
        const json& in = section(root, "input");
        checkKeys(in, "input", {"format", "column", "root_object"});
        const std::string format = in.value("format", "times");
        if (format == "times") cfg.input.format = InputFormat::Times;
        else if (format == "histogram") cfg.input.format = InputFormat::Histogram;
        else if (format == "root") cfg.input.format = InputFormat::Root;
        else fail("input.format must be times / histogram / root");
        cfg.input.column = in.value("column", 0);
        cfg.input.rootObject = in.value("root_object", "");
        if (cfg.input.column < 0) fail("input.column must be >= 0");
        if (cfg.input.format == InputFormat::Root && cfg.input.rootObject.empty())
            fail("input.root_object is required for ROOT input");

        // --- histogram
        const json& h = section(root, "histogram");
        checkKeys(h, "histogram", {"bin_width", "range"});
        if (!h.contains("bin_width")) fail("histogram.bin_width is required");
        cfg.histogram.binWidth = h.at("bin_width").get<double>();
        auto range = readPair<double>(h, "range", {0.0, 0.0});
        cfg.histogram.min = range.first;
        cfg.histogram.max = range.second;
        if (cfg.histogram.binWidth <= 0.0) fail("histogram.bin_width must be > 0");
        if (cfg.input.format == InputFormat::Times && cfg.histogram.max <= cfg.histogram.min)
            fail("histogram.range [min, max] with max > min is required for times");

        // --- estimation
        const json& e = section(root, "estimation");
        checkKeys(e, "estimation", {"min_counts_first_bin", "min_counts_last_bin", "background",
                                    "fit_start_fraction", "fit_start", "fit_end"});
        auto& est = cfg.estimation;
        est.minCountsFirstBin = e.value("min_counts_first_bin", est.minCountsFirstBin);
        est.minCountsLastBin = e.value("min_counts_last_bin", est.minCountsLastBin);
        est.fitStartFraction = e.value("fit_start_fraction", est.fitStartFraction);
        if (e.contains("fit_start")) {
            est.useFitStart = true;
            est.fitStart = e.at("fit_start").get<double>();
        }
        if (!e.contains("fit_end")) fail("estimation.fit_end is required");
        est.fitEnd = e.at("fit_end").get<double>();
        const json& b = section(e, "background");
        checkKeys(b, "estimation.background", {"side", "bins", "skip_bins", "n_sigma", "range"});
        const std::string side = b.value("side", "right");
        if (side != "right" && side != "left") fail("background.side must be right / left");
        est.background.rightSide = side == "right";
        est.background.nBins = b.value("bins", est.background.nBins);
        est.background.skipBins = b.value("skip_bins", est.background.skipBins);
        est.background.nSigma = b.value("n_sigma", est.background.nSigma);
        if (b.contains("range")) {
            std::tie(est.background.from, est.background.to) = readPair<double>(b, "range", {0.0, 0.0});
            if (est.background.to <= est.background.from) fail("background.range must be [from, to] with from < to");
            est.background.useRange = true;
        }
        if (est.useFitStart && est.fitStart >= est.fitEnd) fail("estimation.fit_start must be < fit_end");
        if (est.background.nBins < 2) fail("background.bins must be >= 2");
        if (est.background.skipBins < 0) fail("background.skip_bins must be >= 0");
        if (est.fitStartFraction <= 0.0 || est.fitStartFraction >= 1.0)
            fail("estimation.fit_start_fraction must be in (0, 1)");

        // --- resolution components
        if (!needResolution && root.contains("resolution"))
            std::cerr << "Config warning: \"resolution\" is ignored (formula model without convolution)\n";
        if (needResolution && (!root.contains("resolution") || !root.at("resolution").is_array() ||
                               root.at("resolution").empty()))
            fail("\"resolution\" must be a non-empty array");
        for (std::size_t i = 0; needResolution && i < root.at("resolution").size(); ++i) {
            const json& r = root.at("resolution")[i];
            const std::string where = "resolution[" + std::to_string(i) + "]";
            checkKeys(r, where, {"sigma", "intensity", "sigma_mode", "intensity_mode", "fixed"});
            ResolutionComponent c;
            if (!r.contains("sigma")) fail(where + ".sigma is required (number or \"auto\")");
            if (r.at("sigma").is_string()) {
                if (r.at("sigma").get<std::string>() != "auto") fail(where + ".sigma must be a number or \"auto\"");
                c.sigmaAuto = true;
            } else {
                c.sigma = r.at("sigma").get<double>();
                if (c.sigma <= 0.0) fail(where + ".sigma must be > 0");
            }
            c.intensity = r.value("intensity", 1.0);
            if (r.value("fixed", false)) { c.sigmaMode = ParMode::Fixed; c.intensityMode = ParMode::Fixed; }
            if (r.contains("sigma_mode")) c.sigmaMode = parseMode(r.at("sigma_mode").get<std::string>(), where);
            if (r.contains("intensity_mode")) c.intensityMode = parseMode(r.at("intensity_mode").get<std::string>(), where);
            if (c.sigmaAuto && c.sigmaMode == ParMode::Fixed)
                std::cerr << "Config warning: " << where << " sigma is \"auto\" and fixed\n";
            cfg.resolution.push_back(c);
        }
        if (needResolution)
            checkIntensities<ResolutionComponent>(cfg.resolution, "resolution",
                                                  [](const ResolutionComponent&) { return false; });

        // --- lifetime components (only for the lifetime model)
        if (isFormula && root.contains("lifetimes"))
            std::cerr << "Config warning: \"lifetimes\" is ignored (formula model)\n";
        if (!isFormula && (!root.contains("lifetimes") || !root.at("lifetimes").is_array() ||
                           root.at("lifetimes").empty()))
            fail("\"lifetimes\" must be a non-empty array");
        int nPs = 0;
        for (std::size_t j = 0; !isFormula && j < root.at("lifetimes").size(); ++j) {
            const json& l = root.at("lifetimes")[j];
            const std::string where = "lifetimes[" + std::to_string(j) + "]";
            checkKeys(l, where, {"tau", "intensity", "type", "tau_mode", "intensity_mode", "width", "width_mode",
                                 "comment"});
            LifetimeComponent c;
            if (!l.contains("tau") || !l.contains("intensity")) fail(where + ": tau and intensity are required");
            c.tau = l.at("tau").get<double>();
            c.intensity = l.at("intensity").get<double>();
            if (c.tau <= 0.0) fail(where + ".tau must be > 0");
            if (l.contains("type")) applyLegacyType(l.at("type").get<std::string>(), c, where);
            if (l.contains("tau_mode")) c.tauMode = parseMode(l.at("tau_mode").get<std::string>(), where);
            if (l.contains("intensity_mode")) {
                if (c.pPs) fail(where + ": the p-Ps ratio is always fixed, intensity_mode not allowed");
                c.intensityMode = parseMode(l.at("intensity_mode").get<std::string>(), where);
            }
            if (l.contains("width")) {
                if (c.pPs) fail(where + ": the p-Ps component cannot have a width");
                c.hasWidth = true;
                c.width = l.at("width").get<double>();
                if (c.width < 0.0) fail(where + ".width (sd of ln tau) must be >= 0");
                if (l.contains("width_mode")) c.widthMode = parseMode(l.at("width_mode").get<std::string>(), where);
            } else if (l.contains("width_mode")) {
                fail(where + ": width_mode without width");
            }
            if (c.pPs) {
                ++nPs;
                if (c.intensity <= 0.0) fail(where + ": p-Ps intensity (ratio to o-Ps) must be > 0");
            }
            cfg.lifetimes.push_back(c);
        }
        if (nPs > 1) fail("only one p-Ps (type \"ps\") component is allowed");
        if (!isFormula)
            checkIntensities<LifetimeComponent>(cfg.lifetimes, "lifetime",
                                                [](const LifetimeComponent& c) { return c.pPs; });

        // --- fit
        const json& f = section(root, "fit");
        checkKeys(f, "fit", {"iterations", "partial_variation_percent", "lifetime_limits", "sigma_limits",
                             "offset_window_sigma", "ops_lifetime_threshold", "first_pass_options", "options",
                             "width_limits", "profile_widths"});
        auto& fit = cfg.fit;
        fit.iterations = f.value("iterations", fit.iterations);
        fit.partialVariation = f.value("partial_variation_percent", 100.0 * fit.partialVariation) / 100.0;
        std::tie(fit.tauMin, fit.tauMax) = readPair<double>(f, "lifetime_limits", {fit.tauMin, fit.tauMax});
        std::tie(fit.sigmaMin, fit.sigmaMax) = readPair<double>(f, "sigma_limits", {fit.sigmaMin, fit.sigmaMax});
        fit.offsetWindow = f.value("offset_window_sigma", fit.offsetWindow);
        fit.oPsThreshold = f.value("ops_lifetime_threshold", fit.oPsThreshold);
        fit.firstPassOptions = f.value("first_pass_options", fit.firstPassOptions);
        fit.options = f.value("options", fit.options);
        std::tie(fit.widthMin, fit.widthMax) = readPair<double>(f, "width_limits", {fit.widthMin, fit.widthMax});
        fit.profileWidths = f.value("profile_widths", fit.profileWidths);
        if (fit.widthMin < 0.0 || fit.widthMax <= fit.widthMin) fail("fit.width_limits must be 0 <= low < high");
        if (fit.iterations < 0) fail("fit.iterations must be >= 0");
        if (fit.partialVariation <= 0.0) fail("fit.partial_variation_percent must be > 0");
        if (fit.tauMin <= 0.0 || fit.tauMax <= fit.tauMin) fail("fit.lifetime_limits must be 0 < low < high");
        if (fit.sigmaMin <= 0.0 || fit.sigmaMax <= fit.sigmaMin) fail("fit.sigma_limits must be 0 < low < high");
        if (fit.offsetWindow <= 0.0) fail("fit.offset_window_sigma must be > 0");

        // --- distribution analysis (optional)
        if (root.contains("distribution")) {
            const json& d = section(root, "distribution");
            checkKeys(d, "distribution", {"enabled", "methods", "grid", "bootstrap", "seed", "regions",
                                          "statistics_thresholds", "comment"});
            auto& dc = cfg.distribution;
            dc.enabled = d.value("enabled", true);
            if (d.contains("methods")) {
                dc.nnlsBagged = dc.maxent = false;
                for (const auto& m : d.at("methods")) {
                    const std::string name = m.get<std::string>();
                    if (name == "nnls_bagged") dc.nnlsBagged = true;
                    else if (name == "maxent") dc.maxent = true;
                    else fail("distribution.methods: unknown method \"" + name + "\" (nnls_bagged, maxent)");
                }
                if (!dc.nnlsBagged && !dc.maxent) fail("distribution.methods is empty");
            }
            if (d.contains("grid")) {
                const json& g = d.at("grid");
                if (!g.is_array() || g.size() != 3) fail("distribution.grid must be [tau_min, tau_max, points]");
                dc.gridMin = g[0].get<double>();
                dc.gridMax = g[1].get<double>();
                dc.gridPoints = g[2].get<int>();
            }
            if (dc.gridMin <= 0.0 || dc.gridMax <= dc.gridMin || dc.gridPoints < 10)
                fail("distribution.grid must be 0 < tau_min < tau_max with >= 10 points");
            dc.bootstrap = d.value("bootstrap", dc.bootstrap);
            dc.seed = d.value("seed", dc.seed);
            if (dc.bootstrap < 2) fail("distribution.bootstrap must be >= 2");
            if (d.contains("regions")) {
                for (const auto& r : d.at("regions")) {
                    if (!r.is_array() || r.size() != 2 || r[1].get<double>() <= r[0].get<double>())
                        fail("distribution.regions must be a list of [tau_low, tau_high]");
                    dc.regions.emplace_back(r[0].get<double>(), r[1].get<double>());
                }
            }
            std::tie(dc.countsLow, dc.countsGood) =
                readPair<double>(d, "statistics_thresholds", {dc.countsLow, dc.countsGood});
        }

        // --- output
        const json& o = section(root, "output");
        checkKeys(o, "output", {"directory", "timestamp"});
        cfg.output.directory = o.value("directory", cfg.output.directory);
        cfg.output.timestamp = o.value("timestamp", cfg.output.timestamp);
    } catch (const json::exception& ex) {
        fail(std::string("wrong value type: ") + ex.what());
    }
    return cfg;
}

// ----------------------------------------------------------------------------- writing
std::string Config::toJson() const
{
    using oj = nlohmann::ordered_json;
    auto pair = [](double a, double b) { return oj::array({a, b}); };
    oj root;
    if (!description.empty()) root["description"] = description;

    oj in;
    in["format"] = input.format == InputFormat::Times ? "times" : input.format == InputFormat::Histogram ? "histogram" : "root";
    if (input.format != InputFormat::Root) in["column"] = input.column;
    else in["root_object"] = input.rootObject;
    root["input"] = in;

    oj h;
    h["bin_width"] = histogram.binWidth;
    h["range"] = pair(histogram.min, histogram.max);
    root["histogram"] = h;

    oj e;
    e["min_counts_first_bin"] = estimation.minCountsFirstBin;
    e["min_counts_last_bin"] = estimation.minCountsLastBin;
    oj b;
    b["side"] = estimation.background.rightSide ? "right" : "left";
    b["bins"] = estimation.background.nBins;
    b["skip_bins"] = estimation.background.skipBins;
    b["n_sigma"] = estimation.background.nSigma;
    if (estimation.background.useRange) b["range"] = pair(estimation.background.from, estimation.background.to);
    e["background"] = b;
    e["fit_start_fraction"] = estimation.fitStartFraction;
    if (estimation.useFitStart) e["fit_start"] = estimation.fitStart;
    e["fit_end"] = estimation.fitEnd;
    root["estimation"] = e;

    const bool formula = model.type == ModelType::Formula;
    if (formula) {
        oj m;
        m["type"] = "formula";
        m["formula"] = model.formula;
        m["convolve"] = model.convolve;
        m["points_per_sigma"] = model.pointsPerSigma;
        if (model.step > 0.0) m["step"] = model.step;
        m["window_sigma"] = model.windowSigma;
        oj pars = oj::object();
        for (const auto& p : model.parameters) {
            oj q;
            if (p.estimate.empty()) q["value"] = p.value;
            else {
                std::ostringstream v;
                v << "estimate:" << p.estimate;
                if (p.estimateFactor != 1.0) v << "*" << std::setprecision(15) << p.estimateFactor;
                q["value"] = v.str();
            }
            if (p.hasLimits) q["limits"] = pair(p.lower, p.upper);
            q["mode"] = toString(p.mode);
            q["stage"] = p.stage;
            pars[p.name] = q;
        }
        m["parameters"] = pars;
        root["model"] = m;
    }

    if (!formula || model.convolve) {
        oj res = oj::array();
        for (const auto& r : resolution) {
            oj q;
            if (r.sigmaAuto) q["sigma"] = "auto";
            else q["sigma"] = r.sigma;
            q["intensity"] = r.intensity;
            q["sigma_mode"] = toString(r.sigmaMode);
            q["intensity_mode"] = toString(r.intensityMode);
            res.push_back(q);
        }
        root["resolution"] = res;
    }

    if (!formula) {
        oj life = oj::array();
        for (const auto& l : lifetimes) {
            oj q;
            q["tau"] = l.tau;
            q["intensity"] = l.intensity;
            if (l.pPs) q["type"] = "ps";
            q["tau_mode"] = toString(l.tauMode);
            if (!l.pPs) q["intensity_mode"] = toString(l.intensityMode);
            if (l.hasWidth) {
                q["width"] = l.width;
                q["width_mode"] = toString(l.widthMode);
            }
            life.push_back(q);
        }
        root["lifetimes"] = life;
    }

    oj f;
    f["iterations"] = fit.iterations;
    f["partial_variation_percent"] = 100.0 * fit.partialVariation;
    f["lifetime_limits"] = pair(fit.tauMin, fit.tauMax);
    f["sigma_limits"] = pair(fit.sigmaMin, fit.sigmaMax);
    f["offset_window_sigma"] = fit.offsetWindow;
    f["ops_lifetime_threshold"] = fit.oPsThreshold;
    f["first_pass_options"] = fit.firstPassOptions;
    f["options"] = fit.options;
    f["width_limits"] = pair(fit.widthMin, fit.widthMax);
    f["profile_widths"] = fit.profileWidths;
    root["fit"] = f;

    if (distribution.enabled) {
        oj d;
        oj methods = oj::array();
        if (distribution.nnlsBagged) methods.push_back("nnls_bagged");
        if (distribution.maxent) methods.push_back("maxent");
        d["methods"] = methods;
        d["grid"] = oj::array({distribution.gridMin, distribution.gridMax, distribution.gridPoints});
        d["bootstrap"] = distribution.bootstrap;
        d["seed"] = distribution.seed;
        if (!distribution.regions.empty()) {
            oj regs = oj::array();
            for (const auto& r : distribution.regions) regs.push_back(pair(r.first, r.second));
            d["regions"] = regs;
        }
        d["statistics_thresholds"] = pair(distribution.countsLow, distribution.countsGood);
        root["distribution"] = d;
    }

    oj o;
    o["directory"] = output.directory;
    o["timestamp"] = output.timestamp;
    root["output"] = o;
    return root.dump(2) + "\n";
}

void Config::save(const std::string& path) const
{
    std::ofstream out(path);
    if (!out) throw std::runtime_error("Cannot write config file " + path);
    out << toJson();
}

} // namespace pals
