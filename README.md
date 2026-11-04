# unitroottools

[![R-CMD-check](https://github.com/maarranz/UnitRoot-Tools-R/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/maarranz/UnitRoot-Tools-R/actions/workflows/R-CMD-check.yaml)

Transparent unit-root testing in R: one interface for classical unit-root and
stationarity tests, plus **deterministic-specification strategies** that record
every step of the decision. Developed for *Applied Time Series Econometrics*
(ATSE) at Universidad Carlos III de Madrid; companion to the Python package
UnitRoot-Tools.

> **Status: teaching version (v0.1 in preparation).** The functions listed
> under [What it does](#what-it-does) are implemented and validated. Interfaces
> may still change before the first CRAN release.

## Installation

The package is installed from GitHub. In R (or RStudio), run once:

```r
install.packages("remotes")
remotes::install_github("maarranz/UnitRoot-Tools-R", build_vignettes = TRUE)
```

`build_vignettes = TRUE` also builds the tutorial described below. It needs
[pandoc](https://pandoc.org), which comes with RStudio. If you do not use
RStudio and do not have pandoc, install without the vignette (it can be read
[on the web](https://maarranz.github.io/UnitRoot-Tools-R/vignette.html)):

```r
remotes::install_github("maarranz/UnitRoot-Tools-R")
```

or, equivalently, with [pak](https://pak.r-lib.org):

```r
install.packages("pak")
pak::pak("maarranz/UnitRoot-Tools-R")
```

To update to the latest version, run the same command again. To install a
fixed release instead, add its tag, e.g.
`remotes::install_github("maarranz/UnitRoot-Tools-R@v0.1.1", build_vignettes = TRUE)`
(releases are listed in [NEWS.md](NEWS.md)). The package needs
R 4.1 or later; its non-base dependencies,
[urca](https://CRAN.R-project.org/package=urca) and
[generics](https://CRAN.R-project.org/package=generics), are installed automatically.

## Documentation

- **On the web** (nothing to build): the [API reference](https://maarranz.github.io/UnitRoot-Tools-R/)
  ([source](docs/unitroottools_Documentation.qmd)) and the
  [vignette](https://maarranz.github.io/UnitRoot-Tools-R/vignette.html), a
  tutorial that goes through all the functions with Enders' U.S. real GDP example.
- **In R**: help on each function, e.g. `?adf` or `?ur_strategy`, and, if you
  installed with `build_vignettes = TRUE`, `vignette("unitroottools")`.

## Quick start

```r
library(unitroottools)

gdp  <- read.csv(system.file("fixtures", "enders_rgdp.csv", package = "unitroottools"))
lgdp <- ts(log(gdp$rgdp), start = c(1947, 1), frequency = 4)

adf(lgdp, trend = "ct")                          # ADF test, lags chosen by MAIC
summary(adf(lgdp, trend = "ct"))                 # sources, phi statistics, lag table
ur_strategy(lgdp, "original", lags = 1)          # Enders' testing sequence
ur_strategy(lgdp, "hamilton", trending = TRUE)   # your judgement, then one test
```

Every result states the hypotheses in words, the deterministic terms, how many
lags were used and how they were chosen, the effective sample, the critical
values and where they come from, and the decision.

## What it does

| Function | Purpose |
|---|---|
| `adf()` | Augmented Dickey-Fuller test: τ statistic and the φ statistics, each with critical values and p-value |
| `select_lags()` | Lag order by AIC, BIC, MAIC or MBIC, all candidates on the same sample, with the full table |
| `ur_strategy()` | Choosing the deterministic terms: `"original"` (Dolado-Jenkinson-Sosvilla-Rivero / Enders), `"conservative"`, `"hamilton"`; full history of every step |
| `pp()` | Phillips-Perron Z_τ test, with its p-value, and Z_α with its own critical values |
| `kpss()` | KPSS test, with stationarity as the null hypothesis; critical values and p-value from the exact limiting distribution |
| `dfgls()` | DF-GLS test (Elliott-Rothenberg-Stock), lags chosen by MAIC as recommended by Perron and Qu (2007), with finite-sample critical values and p-value from the package's own simulation |
| `ng_perron()` | Ng-Perron M tests (MZ_α, MZ_t, MSB, MP_T) on GLS-detrended data, lags by MAIC (Perron and Qu, 2007), own finite-sample critical values and p-value |
| `ers_pt()` | Elliott-Rothenberg-Stock point-optimal test P_T, same conventions |
| `boot_ur()` | Bootstrap ADF and union-of-rejections tests through the bootUR package (optional; p-value only) |
| `tidy()`, `glance()` | Results as one-row-per-statistic or one-row-per-test data frames, easy to combine across series |
| `simulate_df()` | Monte Carlo of the Dickey-Fuller distribution, with a plot against N(0,1) |

Details that matter in practice:

- **Critical values.** τ uses MacKinnon's response surfaces at the effective
  sample size, with MacKinnon's (1996) p-values. The φ statistics, DF-GLS,
  the Ng-Perron and ERS tests and Z_α use this package's own finite-sample
  response surfaces in the sample size and the number of lags, estimated by
  Monte Carlo (code in `data-raw/`); KPSS uses its exact limiting
  distribution. Every result records its source and has a p-value.
- **Lags.** Chosen over 0, ..., kmax with kmax = ⌊12(T/100)^¼⌋ (Schwert), all
  candidates on the same sample; the test regression is then estimated on all
  the observations available for the chosen order (`final_sample = "common"`
  keeps the selection sample).
- **Strategies.** The history separates rejections of the unit root from
  rejections that only decide the path (the φ tests and the tests of the
  deterministic terms). No p-value is reported for a whole sequence, because
  its size is not known. Hamilton's rule has no default: you must declare
  whether the series is trending.
- **Data.** Numeric vectors, `ts`, `zoo`, `xts`, data frames, tibbles and
  tsibbles; with several columns, choose one with `var = "name"`. Missing
  values are never dropped silently.
- **Validation.** Statistics are checked against independent `lm()`
  computations and urca, against 120-digit reference calculations shared with
  the Python package, critical values against statsmodels, and the
  strategies against Enders' published GDP example.
- **Long-run variance** (for `pp()` and `kpss()`): Bartlett, Parzen or
  Quadratic Spectral kernel; Schwert's rules, a fixed number of
  autocovariances, or automatic bandwidths (Andrews, Newey-West).
- **Current limits.** The `dfgls()` statistic is computed by urca (checked
  against an independent implementation); everything else is internal. KPSS
  critical values are asymptotic (exact).

## Forthcoming

| Phase | Content |
|---|---|
| v0.2: inference tools (remaining) | Documentation of the new tests; release |
| v0.3: bootstrap | Own bootstrap engine (sieve, wild, block and dependent wild schemes), able to bootstrap any statistic, including a whole testing strategy |
| v0.4: specification uncertainty | Monte Carlo tools (size, power, branch frequencies); union-of-rejections test; bootstrap-calibrated strategies; Hamilton's Case 3 (drift under the null), introduced with the Monte Carlo tools |
| v0.5: robustness | Unit-root tests robust to additive outliers; tests allowing for structural breaks |

## Specification

The econometric conventions (statistics, lag selection, critical values, the
result format and the strategies) are defined in a language-neutral
specification, [`spec/specification.qmd`](spec/specification.qmd), shared with
the Python package UnitRoot-Tools. The two packages are independent
implementations of that specification and are validated against shared
numerical fixtures ([`inst/fixtures`](inst/fixtures)).

## Licence

MIT © Miguel A. Arranz
