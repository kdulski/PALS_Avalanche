# PALS_avalanche

Discrete fit of positron annihilation lifetime (PAL) spectra with ROOT.

## Build

Requirements: ROOT ≥ 6.22, a C++17 compiler. (Boost is no longer needed; the
JSON parser `nlohmann/json` 3.12.0 is bundled in `third_party/`.)

```sh
make                 # -> ./PALS_avalanche
make test            # ROOT-independent tests + benchmark (uses the example data)
# or
cmake -B build && cmake --build build
```

## Usage

### Graphical interface

```sh
./PALS_avalanche_gui [data file] [config.json]      # build: make gui, or CMake (PALS_BUILD_GUI=ON)
```

* left: the spectrum (log y) with the fit, components, fit range (blue) and background window
  (orange), residuals below; mouse modes: **zoom** (ROOT default), **fit range** and
  **background window** (drag over the spectrum); the cursor position and counts are shown below
* right: load data / config, save config; data settings (format, column, ROOT object, bin
  width, range – changes re-read the data); fit start (explicit or fraction of the peak), fit
  end, background window, current estimates; resolution and lifetime components (tau,
  intensity, modes, p-Ps, width), add / remove; fit settings; **Fit** and **Save results**;
  stage progress and the report
* everything set in the GUI is the same `Config` as a file: **Save config** writes a JSON that
  the command-line program reads; settings are validated like a config file
* formula models are shown and fitted, but their formula is edited in the config file
* look: colours, fonts, button shape and where the PAv logo appears (Fit panel, below the
  residuals, watermark behind the plot) are set in `GuiTheme` (`src/gui/PalsGui.h`); the logo
  images are read from `resources/` (source tree, next to the executable, the working
  directory, or `$PALS_RESOURCES`) – without them the GUI simply shows no logo

### Command line

```sh
./PALS_avalanche <data file> [config.json]      # default config: FitDetails.json
./PALS_avalanche ColonCancer_AB_20122018.dat_TDiff.dat config/tumor.json
```

Output (`output.directory`, default `Results/`):

* `<data stem>_<timestamp>.txt` – estimates, fit stages, results with errors, all parameters
* `<data stem>_<timestamp>.root` – `data` (histogram with the fit attached), `fit`,
  `component_<j>`, `background`, `residuals`, `config` (the JSON used)

Old `FitDetails` files can be converted with
`tools/fitdetails_to_json.py FitDetails_old new.json`.

## Configuration

See `config/tumor.json`. JSON may contain `//` comments. Unknown keys are reported.

| key | default | meaning |
|---|---|---|
| `input.format` | `times` | `times` (list of time differences), `histogram` (text, one count per bin), `root` |
| `input.column` | 0 | column of a text file (0-based; oscilloscope files: 1) |
| `input.root_object` | – | `"dir/histogram"` in a ROOT file |
| `histogram.bin_width` | required | ns; for ROOT input: an integer multiple of the original bin width (rebinning) |
| `histogram.range` | required for `times` | `[low edge, high edge]` ns; for `histogram` only the low edge is used |
| `estimation.min_counts_first_bin` / `min_counts_last_bin` | 8 | spectrum = first/last bin above these counts |
| `estimation.background.side` | `right` | which end of the spectrum gives the background |
| `estimation.background.bins` / `skip_bins` | 100 / 30 | averaged bins / bins skipped at the end |
| `estimation.background.n_sigma` | 3 | background limits: estimate ± n_sigma × its uncertainty |
| `estimation.fit_start_fraction` | 0.3 | fit starts where (N − B) on the leading edge drops below this fraction of the peak |
| `estimation.fit_start` | – | explicit start of the fit [ns] (instead of `fit_start_fraction`) |
| `estimation.background.range` | – | explicit background window `[from, to]` ns (instead of side / bins / skip) |
| `estimation.fit_end` | required | ns |
| `resolution[]` | required | `{"sigma": 0.1 or "auto", "intensity": 1, "sigma_mode": "free", "intensity_mode": "free"}`; `"fixed": true` = old "fix the resolution" |
| `lifetimes[]` | required | `{"tau": 0.4, "intensity": 0.3, "type": "nf"}` or with `tau_mode` / `intensity_mode` |
| `fit.iterations` | 3 | iterations after the first stage |
| `fit.partial_variation_percent` | 5 | window of "partial" parameters |
| `fit.lifetime_limits` | [0.08, 142] | limits of free lifetimes, ns |
| `fit.sigma_limits` | [0.01, 5] | limits of free sigmas, ns |
| `fit.offset_window_sigma` | 5 | t0 ∈ peak ± k·sigma |
| `fit.ops_lifetime_threshold` | 0.7 | components above it are o-Ps (only with a p-Ps component) |
| `fit.first_pass_options` / `fit.options` | `RMW` / `RM` | ROOT fit options of stage 0 / later stages (e.g. `RML` for Poisson likelihood) |
| `output.directory` / `output.timestamp` | `Results` / true | |

### Model given by a formula

```json
"model": {
  "type": "formula",
  "formula": "[A]*([I1]*exp(-x/[tau1])/[tau1] + (1-[I1])*exp(-x/[tau2])/[tau2])",
  "convolve": true,
  "parameters": {
    "A":    {"value": "estimate:area"},
    "I1":   {"value": 0.7, "limits": [0, 1]},
    "tau1": {"value": 0.2, "limits": [0.08, 142], "stage": 1},
    "tau2": {"value": 2.0, "limits": [0.08, 142], "mode": "partial", "stage": 1}
  }
}
```

* `formula` – TFormula expression in `x` with named parameters `[name]` (numbered `[0]`
  are called `p0`, ...). It is a **density in counts per ns**; the program multiplies it by
  the bin width and adds the background (fitted as before).
* `convolve: true` – `x` is the time since t0 (the formula is taken as 0 for x < 0) and the
  formula is convolved numerically with the resolution function from `"resolution"`
  (sigma, t0, fractions fitted as in the lifetime model). Grid step
  `min(sigma)/points_per_sigma` (default 5) or `step` [ns]; Gaussian truncated at
  `window_sigma` (default 8). Quadratic interpolation on the grid, exact Gaussian moments:
  relative accuracy ~1e-5 for tau = 0.125 ns, sigma = 0.12 ns.
* `convolve: false` – `x` is the time itself and the formula must contain the resolution.
  The function `pals_emg(x, sigma, t0, tau)` (exponential convolved with a Gaussian, unit
  area) can be used for that. `"resolution"` is then ignored.
* parameters: `value` – number or `"estimate:<name>"` / `"estimate:<name>*factor"` with
  `background, background_error, area, peak_position, peak_counts, sigma, bin_width,
  fit_start, fit_end`; `limits` (optional); `mode` free / fixed / partial; `stage` – first
  fit stage in which the parameter is released (default 0).

Examples: `config/tumor_formula.json` (same model as `tumor.json`, numerical convolution),
`config/tumor_formula_emg.json` (no convolution, `pals_emg` in the formula).

### Width of a component (log-normal distribution of lifetimes)

```json
{"tau": 2.0, "intensity": 0.1, "type": "nf", "width": 0.05, "width_mode": "free"}
```

`width` = standard deviation of ln(tau) of a log-normal distribution whose mean is `tau`
(sd of tau ≈ tau·sqrt(exp(width²) − 1)); integrated with 32-point Gauss–Hermite quadrature.
Released together with the lifetimes (stage 1), limits `fit.width_limits` (default [0, 1]);
a free width starts at ≥ 0.02 (at 0 the model does not depend on it to first order).
With `fit.profile_widths` (default true) the report gives a 68 % profile-likelihood interval
(the fit statistic + 1, the width fixed on a scan and all other parameters refitted) — the
Hessian error is not reliable near 0.

In simulations (resolution FWHM 0.27 ns, o-Ps 20 %, 2.0 + 2.5 ns or a log-normal with sd 0.25 ns)
the width was 0.105 ± 0.019 (truth 0.111) at 10⁷ counts and undetermined at 10⁶; a single
component gave 0.02 ± 0.03. The report adds a note on this whenever widths are used.

### Lifetime distribution (optional, lifetime model only)

```json
"distribution": {
  "methods": ["nnls_bagged", "maxent"],
  "grid": [0.08, 15.0, 160],
  "bootstrap": 30,
  "seed": 12345,
  "regions": [[1.0, 6.0]],
  "statistics_thresholds": [5e6, 1e7]
}
```

After the discrete fit, with its resolution, t0, background, fit range and Pearson weights:

* **NNLS bagged** – non-negative least squares (Lawson–Hanson, Householder QR) on a log-spaced
  grid of lifetimes, averaged over `bootstrap` Poisson replicas of the NNLS-fitted spectrum;
  spreads (intensity, centroid, 16–84 % of the width, fraction of replicas with ≥ 2 maxima).
* **MaxEnt** – χ²/2 − α·S with a flat default model (damped Newton), α such that
  χ² = χ²(NNLS) + sqrt(2·bins).
* **Comparison** for each region: intensity, centroid, width and number of maxima from the
  sparse (NNLS) and the smooth (MaxEnt) prior; what differs is flagged as depending on the prior.
  Agreement is necessary but not sufficient (two close components and a continuous
  distribution with the same mean and width give the same result).
* The report starts with a note on statistics: in the simulations these methods separated a
  broadened o-Ps from one exponential reliably only from about 5·10⁶ – 10⁷ counts
  (`statistics_thresholds`), and the spectrum's counts are compared with that.
* ROOT file: `distribution_nnls_bagged` (+ `_q16`, `_q84` bootstrap band) and
  `distribution_maxent` – density per ln τ, unit total intensity, log-spaced bins in τ.

Example: `config/tumor_distribution.json` (ColonCancer, ~7·10⁵ counts – below the threshold,
which the report states).

### Component modes

Every lifetime and intensity is `free`, `fixed` or `partial` (± `partial_variation`
around the current value, re-centred in each iteration). Intensities in mode
`free` form one group described by the hyperspherical parametrisation and share
what is left after the `fixed`/`partial` ones; those are absolute intensities.

| old type | tau | intensity |
|---|---|---|
| `nf` | free | free |
| `f` | fixed | fixed |
| `pf` | partial | free |
| `lpf` | partial | partial |
| `lf` | fixed | free |
| `lff` | fixed | partial |
| `ps` | partial | ratio I(p-Ps)/I(o-Ps), always fixed |

Note: in the old code the free intensity of `pf` / `lf` was an absolute
intensity in [0, 1] outside the normalised group; now it belongs to the group,
so the intensities always sum to 1.

### Fit procedure

* stage 0 (`first_pass_options`): all lifetimes fixed; background, area, resolution
  and the free intensities fitted;
* stages 1..`iterations` (`options`): lifetimes and intensities released according
  to their modes.

In general every parameter has a mode and a release stage; `partial` parameters are
never released in stage 0. The report notes parameters that ended at a limit. Status
4000 of the fit only means that IMPROVE (option `M`) found no better minimum.

Initial values: background = mean of the background window, area = background-subtracted
counts of the spectrum, t0 = peak position, sigma = from config or (`"auto"`) from the
left half-width of the peak, lifetimes/intensities = from config.

## Model

N(t) = B + A·Δt·Σᵢ wᵢ Σⱼ Iⱼ f(t; σᵢ, t0ᵢ, τⱼ), where f is an exponential convolved with a
Gaussian (unit area), Δt the bin width and A the number of counts. Resolution weights and
free intensities: hyperspherical parametrisation (`src/Intensities.h`) with one parameter
per component (r, φ₁ … φₙ₋₂, last one unused and fixed). Intensities and their errors in the
report are computed from the full covariance matrix.

## Code

| file | ROOT | content |
|---|---|---|
| `src/Config.*` | – | JSON config |
| `src/Spectrum.*` | – | reading text data, binning, preliminary estimates |
| `src/Intensities.h` | – | hyperspherical parametrisation |
| `src/LifetimeModel.*` | – | lifetime model (functor usable by TF1) |
| `src/FormulaModel.*`, `src/Convolution.*` | – | formula model, numerical convolution |
| `src/Emg.h` | – | `pals_emg`, also declared in the ROOT interpreter |
| `src/FitParameters.h`, `src/ModelSetup.*` | – | parameter specs (initial values, limits, stages) |
| `src/RootFormula.*` | yes | TFormula compilation |
| `src/Quadrature.h` | – | Gauss–Hermite quadrature (component widths) |
| `src/Nnls.*` | – | Lawson–Hanson NNLS, Cholesky |
| `src/Distribution.*` | – | NNLS bagged, MaxEnt, comparison, report section |
| `src/TextReport.*` | – | txt report, error propagation |
| `src/Fitter.*` | yes | TF1/TH1D, staged fit of any model |
| `src/FitSession.*` | yes | one analysis (config, data, estimates, fit, report, files) – CLI and GUI |
| `src/gui/*` | yes (Gui) | `PALS_avalanche_gui` |
| `resources/` | – | PAv logo: `pav_logo.svg` (source), PNGs for the GUI and the window icon |
| `src/RootIO.*` | yes | ROOT input / output |
| `tests/` | – | tests against the original model functions, benchmark |
