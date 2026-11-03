# Shared numerical fixtures

JSON fixtures shared with the Python package UnitRoot-Tools. Both
implementations must reproduce every fixture within its stated tolerance
(see the Validation standard in `spec/specification.qmd`).

Each file contains:

- `id`, `description`
- `source`: provenance of the data and of the expected values, with a hash of the data source
- `data`: the series
- `settings`: test or strategy options
- `expected`: statistics, critical values, decisions, history
- `tolerance`: `rtol` and `atol`, with a justification when looser than the default

## Files

- `enders_rgdp.csv`, `enders_rgdp.json`: Enders' U.S. real GDP example, the
  published fixture for the specification strategies.
- `python/`: high-precision ADF references copied from the Python package
  (see `python/README.md`).
- `dfgls_response_surfaces.json`: the DF-GLS critical-value and p-value
  surfaces (Decision D11), for use by any implementation. Fields: `model` and
  `usage` (how to evaluate critical values and p-values, and what to do
  outside the simulated design), `statistic` (exact definition), `source`
  (simulation and validation, code in `data-raw/`), `trends$c` and
  `trends$ct` (`probs`, `k`, `design`, `intercept`, and `slopes[k][p]` =
  coefficients of 1/T, 1/T², 1/T³), and `reference_values` (critical values
  and p-values computed by R at 66 points, which an implementation reading
  the file should reproduce to 1e-12). Numbers are written with 17
  significant digits, so they reproduce R's coefficients exactly.
