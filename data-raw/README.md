# data-raw: DF-GLS response surfaces (Decision D11)

Reproduces `R/sysdata.rda` (object `dfgls_rs`). Run from this folder; needs Rcpp.

1. `dfgls_simulate.R` — DF-GLS statistic under the null (random walk, N(0,1)
   innovations) for T = 20, ..., 1000 and k = 0, ..., min(12, Schwert's kmax);
   trends c and ct; 40 experiments x 10,000 replications per (T, trend), all
   lag orders on the same paths; quantiles at 117 probabilities. Output
   `dfgls_quantiles.csv.gz`. About 14 minutes on 2 cores.
2. `dfgls_simulate_ext.R` — the same for k = 13, ..., 21 at T = 150, ..., 1000.
   Output `dfgls_quantiles_ext.csv.gz`.
3. `Rscript dfgls_response_surfaces.R ../R/sysdata.rda` — weighted least
   squares, for each trend and probability: common asymptotic value plus a
   cubic in 1/T for each k, on T >= 25.
4. `dfgls_validate.R` — size of the 1/5/10% critical values in fresh
   simulations at (T, k) not used in the fit.

`dfgls_kernel.cpp` computes the statistic exactly as `dfgls()` does (checked
to 4e-15 against `dfgls()` and an independent implementation). Seeds are
fixed (L'Ecuyer-CMRG streams), so the files are reproducible.
