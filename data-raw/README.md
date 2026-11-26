# data-raw: calibration of critical values and p-values

Reproduces `R/sysdata.rda`. Run everything from this folder; the simulations need Rcpp. `make_sysdata.R` assembles `R/sysdata.rda` from the `.rds` files the scripts below produce.

## DF-GLS response surfaces (Decision D11): `dfgls_rs`

1. `dfgls_simulate.R` — DF-GLS statistic under the null (random walk, N(0,1)
   innovations) for T = 20, ..., 1000 and k = 0, ..., min(12, Schwert's kmax);
   trends c and ct; 40 experiments x 10,000 replications per (T, trend), all
   lag orders on the same paths; quantiles at 117 probabilities. Output
   `dfgls_quantiles.csv.gz`. About 14 minutes on 2 cores.
2. `dfgls_simulate_ext.R` — the same for k = 13, ..., 21 at T = 150, ..., 1000.
   Output `dfgls_quantiles_ext.csv.gz`.
3. `Rscript dfgls_response_surfaces.R dfgls_rs.rds` — weighted least
   squares, for each trend and probability: common asymptotic value plus a
   cubic in 1/T for each k, on T >= 25.
4. `dfgls_validate.R` — size of the 1/5/10% critical values in fresh
   simulations at (T, k) not used in the fit.

`dfgls_kernel.cpp` computes the statistic exactly as `dfgls()` does (checked
to 4e-15 against `dfgls()` and an independent implementation). Seeds are
fixed (L'Ecuyer-CMRG streams), so the files are reproducible.

## KPSS limiting distributions (Decision D14): `kpss_limit`

- `kpss_limit.R` — eigenvalues of the covariance kernels of the first- and
  second-level Brownian bridges (Nystrom method, 2000 points). Output
  `kpss_limit.rds`. Probabilities are computed from them in the package by
  Imhof's formula (`kpss_upper_prob()` in `R/kpss.R`).
- `kpss_simulate.R` (kernel `kpss_kernel.cpp`) — independent check by
  simulation (T = 2000, 1,000,000 replications); agrees with the exact
  quantiles to simulation error.

## φ statistics and the normalized bias (Decision D16): `phi_rs`, `rho_rs`

1. `phi_simulate.R` (kernel `phi_kernel.cpp`) — φ1, φ2, φ3 (and τ2, τ3 as a
   check) of the ADF regression under the null, same design as DF-GLS
   (T = 20, ..., 1000; k = 0, ..., min(12, kmax), plus 13, ..., 21 for
   T >= 150); 40 x 10,000 replications. Output `phi_quantiles.csv.gz`.
   About 50 minutes on 2 cores.
2. `rho_simulate.R` (kernel `rho_kernel.cpp`) — the Dickey-Fuller normalized
   bias T(rho - 1) without lags (k = 0), used for the critical values of
   Phillips-Perron's Z_alpha; 40 x 25,000 replications. Output
   `rho_quantiles.csv.gz`.
3. `Rscript new_response_surfaces.R` — fits both with `fit_surface()`
   (`surface_fit.R`, as for DF-GLS) and writes `phi_rs.rds`, `rho_rs.rds`
   (and `gls_rs.rds`, below).

## Ng-Perron and ERS statistics (Decision D15): `gls_rs`

1. `gls_simulate.R` (kernel `gls_kernel.cpp`) — MZa, MZt, MSB, MPT and PT
   under the null, trends c and ct, same (T, k) design; 20 x 10,000
   replications. Output `gls_quantiles.csv.gz`.
2. `Rscript new_response_surfaces.R` — these are *tabulated*, not fitted:
   their tails explode at small T with many lags (the AR spectral estimate),
   so `table_calibration()` stores the simulated quantiles by (k, T) and the
   package interpolates them in 1/T, with a limit fitted on the large-T cells.
   Writes `gls_rs.rds`. The out-of-sample size check of the 0.2.0 release
   (5% size between 4.0 and 5.8%) was run interactively; its script was not
   kept.

## Union of rejections (0.3.4, D22): `union_rs`

- `union_simulate.R` (kernel `union_null_stats()` in `tau_kernel.cpp`) — the
  four statistics (ADF-OLS and DF-GLS, with constant and with trend) on the
  same random-walk paths, normalised by the package's own 5% critical values
  (`tau_rs`, `dfgls_rs`); quantiles of their minimum U. Same (T, k) design as
  τ; 25 x 10,000 replications per T, seed 20261126 (57 minutes on 2 cores).
  The kernel matches `adf()` and `dfgls()` to 2e-13.
- `union_fit.R` → `union_rs.rds` (largest residual at 1/5/10% 0.042, RMS
  ≤ 0.0044). **Depends on `tau_rs` and `dfgls_rs`: rerun the simulation if
  either is refitted.**
- `union_validate.R` — out-of-sample size (4.7%-5.3% at 5%, 14 cells,
  20,000 replications) and power against each component and `ur_strategy()`.
- `union_export_json.R` → `inst/fixtures/union_response_surface.json`.

## ADF τ statistics (0.3.3, D21): `tau_rs`

- `tau_simulate.R` (kernel `tau_kernel.cpp`) — τ1, τ2, τ3 under the null,
  same (T, k) design as φ; 40 x 10,000 replications per T, seed 20261125
  (34 minutes on 2 cores). Output `tau_quantiles.csv.gz`. The kernel also
  computes Engle-Granger residual-based statistics (`nvars > 1`), for
  cointegrationtools; it was checked against `adf()` (1e-14) and against
  `lm()` followed by an ADF regression (exact).
- Fitted in `new_response_surfaces.R` (τ section) → `tau_rs.rds`. The τ2/τ3
  quantiles of `phi_quantiles.csv.gz` (independent seeds) agree with this run
  (sd of the standardised differences 1.005).
- `tau_validate.R` — out-of-sample size check on fresh seeds, at T = 33,
  60, 137, 290, 610 (between grid points) with up to k = 20, 20,000
  replications per cell, using both the critical values and p ≤ level:
  5% size between 4.72% and 5.42% (simulation s.e. 0.15%; 39 cells), the two
  rules agreeing exactly. It also compares the k = 0 values with MacKinnon
  (2010): differences at most 0.008 (ct, N = 24), otherwise below 0.005.

## Assembly

`Rscript make_sysdata.R` writes `../R/sysdata.rda` with every object above.

## Platform and requirements

- Regeneration only: R with Rcpp and a C++ compiler (Xcode command-line
  tools on macOS, Rtools on Windows). Installing or using the package needs
  neither.
- The simulation scripts parallelise with `parallel::mclapply()`, which forks
  processes; on Windows it runs on one core only (set `mc.cores = 1`, or run
  on macOS or Linux). Each job has its own L'Ecuyer-CMRG stream, so results
  do not depend on the number of cores.
- The files in this folder were produced on Linux (R 4.3, 2 cores).
