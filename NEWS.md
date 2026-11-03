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
