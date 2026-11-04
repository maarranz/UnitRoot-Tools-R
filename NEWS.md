# unitroottools (development version)

Fixes from the audit of 0.2.0 (4 November 2026).

* **p-values from the package's own calibrations** (DF-GLS, φ, Z_alpha,
  Ng-Perron, ERS) could be wrong far in the tail: the local cubic used to
  interpolate between simulated quantiles was extrapolated beyond them and
  could turn back, giving p = 1 for clearly stationary data that the test
  rejected (e.g. Ng-Perron MZ_t at T = 250). The distribution function is
  now interpolated linearly in qnorm(p) between the 117 simulated quantiles
  (monotone, exact at them). Beyond the simulated range the p-value is
  reported as the bound 0.001 or 0.999, printed "< 0.001" or "> 0.999",
  with a note in the result. Inside the range, p-values change only in the
  third or fourth decimal. The shared DF-GLS JSON file and its Python reader
  use the same rule.
* `kpss()` no longer fails for very large statistics (e.g. a trending series
  with `bandwidth = 0`, where the numerical integration did not converge):
  beyond the 0.001 and 0.999 quantiles of the limit the p-value is reported
  as the bound, as for the other calibrations. A constant series (or an
  exact linear trend with `trend = "ct"`) now stops with a clear error,
  instead of a statistic made of rounding noise or an integration error.
  The source is described as the limiting distribution "computed
  numerically" rather than "exact".

# unitroottools 0.2.0

* The φ statistics of `adf()` and `ur_strategy()` now take their critical
  values and p-values from the package's own finite-sample response
  surfaces in the sample size and the number of lags (decision D16),
  replacing the Dickey-Fuller (1981) tables. Without lags they agree with
  the tables (e.g. φ3 at T = 100: 8.76/6.48/5.47 against 8.73/6.49/5.47);
  with lags they also account for k. φ statistics now have p-values
  (shown by `summary()`, `tidy()` and the strategy history).
* `pp()`'s Z_alpha now has critical values and a p-value, from the
  finite-sample distribution of the Dickey-Fuller normalized bias (D16).
* New `ng_perron()` (MZ_alpha, MZ_t, MSB, MP_T) and `ers_pt()` (ERS
  point-optimal P_T), on GLS-detrended data with lags chosen as for
  `dfgls()` (decision D15). Critical values and p-values come from the
  package's own simulation over T and k: the simulated quantiles are
  tabulated and interpolated in 1/T, because at small T with many lags
  their tails explode (e.g. the 5% quantile of MZ_alpha at T = 25 with 6
  lags is about -800) and no smooth surface follows them. Out of sample, the
  5% tests reject between 4.0% and 5.8% of the time under the null. With a
  trend, the limits differ from Ng and Perron's (2001) Table 1 (5% MZ_alpha
  -16.6 against -17.3, MZ_t -2.85 against -2.91); an independent simulation
  and urca's DF-GLS agree with the package.
* When the lag order is large for the sample size, `ng_perron()` and
  `ers_pt()` warn that the null distribution is very heavy-tailed (5%
  critical value of MZ_alpha more than twice its no-lag value) and the test
  can hardly reject; `adf()` warns that its φ critical values are accurate
  only to about 5% at T <= 30 with 7 or more lags.
* New `boot_ur()`: bootstrap ADF and union-of-rejections tests through the
  bootUR package (in Suggests; decision D17). The result is a `ur_test`
  whose decision uses the bootstrap p-value.
* New long-run-variance module (decision D13) for `pp()` and `kpss()`:
  Bartlett, Parzen and Quadratic Spectral kernels (`kernel =`), and
  automatic bandwidths `bandwidth = "andrews"` (Andrews, 1991) and
  `"newey-west"` (Newey and West, 1994), besides Schwert's rules and fixed
  values. Defaults are unchanged.
* `kpss()` is computed internally, no longer by urca, and now reports a
  p-value. Its critical values and p-value come from the exact limiting
  distribution (decision D14); Kwiatkowski et al.'s simulated table, used
  before, differs in the third decimal (e.g. 0.146 instead of 0.148 at 5%
  with a trend).
* `pp()` also reports the normalized-bias statistic Z_alpha (in
  `other_statistics`).
* Results of `pp()` and `kpss()` describe the long-run variance used
  ("Bartlett kernel, 4 autocovariances (short rule, ...)").
* The DF-GLS response surfaces are also available as a language-neutral
  file, `inst/fixtures/dfgls_response_surfaces.json`, with reference values,
  so that the Python package can share the same calibration.

# unitroottools 0.1.1

## Bug fix

* `adf()`, `pp()`, `kpss()` and `dfgls()` returned their result invisibly
  since 0.1.0, so calling them at the console (or in the vignette) printed
  nothing. Results are returned visibly again; a test guards against this.

## Changes from the Phase 2A reconciliation with the Python implementation

These follow decisions D1, D3 and D7-D10 of the specification
(`spec/specification.qmd`).

* After lag selection, `adf()` now re-estimates the test regression on all
  the observations available for the chosen order (`final_sample =
  "maximal"`, the new default), as the Python package, ARCH, bootUR and
  EViews do. `final_sample = "common"` keeps the previous behaviour. Results
  with selected lags can therefore change slightly; results with fixed
  `lags` do not.
* `ur_strategy()` follows the same rule (all its regressions share one
  sample) and gains a `final_sample` argument.
* The default `kmax` is capped for very short series (T <= 21) so that every
  candidate regression keeps two residual degrees of freedom; a note reports
  the cap. A `max_lags` above the feasible maximum is an error.
* `adf()` reports MacKinnon's (1996) finite-sample p-value
  (via `urca::punitroot()`); the decision is still taken from the critical
  value, and the two agree to four decimals at the tabulated levels.
* `dfgls()` selects the number of lags automatically when `lags` is not
  given: MAIC (or another criterion) on the OLS-detrended series, as
  recommended by Perron and Qu (2007).
* New `tidy()` and `glance()` methods for tests and strategies (generics
  from the lightweight `generics` package; broom is not needed).
* New tests against 120-digit reference calculations shared with the Python
  package (`inst/fixtures/python/`).

## Other changes

* The history of `ur_strategy()` now reports a p-value for each step
  (MacKinnon, 1996, for the tau tests; Student t for the conditional and
  auxiliary t tests; none for the phi tables).
* `plot()` of `simulate_df()` results is readable at small sizes, and the
  vignette draws its figures larger.
* The vignette covers p-values, automatic lag selection in `dfgls()`, the
  per-step p-values of strategies, and `tidy()`/`glance()`.

## Phase 2B audit

* The statistics of `kpss()` and `dfgls()` (computed by urca) are now tested
  against independent implementations of the published formulas.
* Found two problems in the urca backend, both fixed below: `dfgls()`
  critical values were too liberal in finite samples, and urca's
  Phillips-Perron Z_tau departs from the published formula in small samples.
* `dfgls()` critical values and p-values now come from the package's own
  finite-sample response surfaces (decision D11), estimated from a Monte
  Carlo simulation over sample size and lag order (code in `data-raw/`).
  urca's critical values, used before, rejected a true unit root 8% of the
  time at the 5% level with T = 100 and a constant; the new ones give
  4.8-5.3% in out-of-sample checks.
* `pp()` is now computed internally (decision D12), following Phillips and
  Perron (1988) and Hamilton (1994, eq. 17.6.8); urca's `ur.pp()` is no
  longer used. Values differ from earlier versions by O(1/T) (about 0.1 at
  T = 30, 0.02 at T = 264). `pp()` now reports a p-value (MacKinnon, 1996)
  and the long-run variance it used.

# unitroottools 0.1.0

## Robustness fixes after the Phase 1 audit

* Inputs with more than two dimensions, and matrix-valued data-frame columns,
  are rejected instead of being flattened into one apparent series.
* `max_lags`, `lags` and fixed `bandwidth` values must be single non-negative
  whole numbers; values such as 2.9 or -0.5 are errors instead of being
  silently truncated.
* `level` is checked when a test is called, with a message listing the
  available levels (0.01, 0.05, 0.10; KPSS also 0.025).
* `ur_strategy()` checks that `allow_no_constant` is `TRUE` or `FALSE`.
* Every regression needs at least two residual degrees of freedom, as the
  specification requires.
* `validate_ur_test()` now enforces the whole result contract, including the
  stored decision, canonical trend codes, integer fields, finite accompanying
  statistics and their critical values. The constructor no longer truncates
  non-whole `lags`, `nobs` or `sample`.
* `print()` on a test result now shows the source of the critical values.
* Maintainer email set; specification and README roadmap made consistent
  (Hamilton's Case 3 in v0.4; `tidy()`/`glance()` marked as planned).
