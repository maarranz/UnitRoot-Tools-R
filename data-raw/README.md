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

## Assembly

`Rscript make_sysdata.R` writes `../R/sysdata.rda` with every object above.
