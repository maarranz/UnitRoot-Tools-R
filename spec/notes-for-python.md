# Notes for the Python package (UnitRoot-Tools), from the R work

3 November 2026. Findings from Phase 2A (reconciliation with the specification) and Phase 2B (econometric audit) of UnitRoot-Tools-R that concern the Python implementation. The specification (`specification.qmd`) lists them as known divergences of the Python implementation. Ordered by importance.

## 1. DF-GLS critical values are conservative when lags are used (ARCH 8.0.0)

ARCH's DF-GLS response surfaces depend only on the effective sample N = T − k − 1. In finite samples the DF-GLS quantiles also depend on k, and in the opposite direction: more lags give a lighter left tail, so the critical value should be *less* negative. Indexing by N makes it *more* negative as k grows, so the test is conservative and loses power.

Simulated size of ARCH's 5% critical value (40,000 replications, random-walk null):

| trend | T | k | size |
|---|---|---|---|
| c | 50 | 0 | 4.8% |
| c | 50 | 4 | 3.4% |
| c | 50 | 8 | 2.4% |
| c | 100 | 10 | 3.7% |
| ct | 50 | 4 | 2.6% |
| ct | 50 | 8 | 1.5% |
| ct | 100 | 10 | 2.8% |
| ct | 250 | 12 | 4.0% |

With k = 0 ARCH is fine (4.6–4.8%). MAIC usually selects k > 0 in applied work, so this matters.

R (Decision D11) uses its own surfaces q_p(T, k) = θ∞ + a_k/T + b_k/T² + c_k/T³ (one cubic per k; T = 25, …, 1000; k = 0, …, 21; 117 probabilities); out-of-sample size 4.8–5.3%. Simulation code and output are in `data-raw/` (CSV of simulated quantiles with standard errors, language-neutral). **Suggestion:** Python reads the same coefficients from `inst/fixtures/dfgls_response_surfaces.json` (R repo), so both implementations share one calibration. The file says how to evaluate critical values and p-values and what to do outside the design, and includes 66 reference values; a 30-line NumPy reader reproduces them to 3 × 10⁻¹⁴.

Found on the way: in small samples the quantiles for k = 2j and k = 2j + 1 are nearly equal (a paired pattern, confirmed with independent code). Polynomials in k cannot fit it; per-k surfaces can.

*Reproduction:* ARCH 8.0.0 `DFGLS(y, trend, lags=k).critical_values["5%"]` gives the critical value Python uses for each (T, k); the null distribution was simulated with `data-raw/dfgls_kernel.cpp` (seed 4242), and the size is the share of simulated statistics below ARCH's value.

## 2. τ p-values: asymptotic (MacKinnon 1994) against finite-sample critical values

Python reports asymptotic p-values with finite-sample (MacKinnon 2010) critical values, so the two disagree near the boundary. A statistic exactly at the 5% critical value gets p = 0.030 (ct, N = 25), 0.036 (c, N = 25), 0.044 (ct, N = 100). Specification D8: MacKinnon (1996) finite-sample p-values at the effective N. These agree with the 2010 critical values to four decimals (p = 0.0100, 0.0500, 0.1000 at the critical values for N = 25, …, 262).

## 3. Default kmax (D3)

- Rounding: Python uses ⌈12(T/100)^¼⌉ (ARCH's convention); Schwert (1989) and the specification use the integer part.
- Cap: Python caps at ⌊(N − 3)/2⌋, which leaves one residual degree of freedom in the largest ADF regression. The specification caps at ⌊(T − 4 − d)/2⌋ (d = number of deterministic terms), so every candidate regression keeps two, and warns when the cap binds (only for T ≤ 21).

## 4. Plain AIC/BIC (D7): admissible, but document it

Python computes AIC/BIC from the detrended regression used by MAIC (so AIC equals MAIC minus the τ term). The specification's default is the ADF regression with its deterministic terms, which is what users expect "AIC" to mean (urca and ARCH's autolag do this). Python's version is recorded as an admissible alternative and should be documented as such.

## 5. Already in agreement (no action)

- τ critical values: R and ARCH identical to 1.8 × 10⁻¹⁵ (n, c, ct; N = 3, …, 1000 and ∞).
- MAIC/MBIC formula and conventions (OLS detrending, common sample, SSR/n, log n, ties to the smaller k).
- Final regression on the maximal sample after selection (R changed to match Python, D1).
- DF-GLS lag selection by Perron–Qu (R changed to match Python, D9).
- φ tables (Enders Table B) and their lookup rule; the Enders strategy logic.
- `validation/enders/scaled_references.json` (120 digits): R matches all 18 cases to < 3 × 10⁻¹⁵.

## 6. Not yet in Python

- Hamilton's rule-based strategy (specification, section "Hamilton (rule-based)").
- Per-step p-values in the strategy history.
- Phillips–Perron: if added through ARCH, note that ARCH follows the published formula (it agrees with R's internal Z_τ to 2 × 10⁻¹³); urca does not.
