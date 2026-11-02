# unitroottools

[![R-CMD-check](https://github.com/maarranz/UnitRoot-Tools-R/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/maarranz/UnitRoot-Tools-R/actions/workflows/R-CMD-check.yaml)

Transparent unit-root testing in R: one interface for classical unit-root and
stationarity tests, plus **deterministic-specification strategies** that record
every step of the decision.

> **Status: early development, not yet usable.** The package structure, the common
> specification and the result format are in place; the tests themselves are being
> implemented. A first teaching release (v0.1) is planned for the *Applied Time
> Series Econometrics* course.

## Why another unit-root package?

Unit-root tests are well covered in R (notably by
[urca](https://CRAN.R-project.org/package=urca) and
[bootUR](https://CRAN.R-project.org/package=bootUR)). `unitroottools` builds on
them and adds:

- **One result format for every test**, whatever computes it: the hypotheses in
  words, the deterministic specification, lags and how they were chosen, the
  effective sample, critical values and their source, the decision, and the
  engine used.
- **Specification strategies** for choosing the deterministic terms: the
  Dolado–Jenkinson–Sosvilla-Rivero / Enders sequence, a conservative variant, and
  Hamilton's rule-based choice, each with a full history of tests and branch
  decisions.
- **Transparent lag selection** (AIC, BIC, MAIC, MBIC on a common sample), showing
  the whole criterion table, not only the selected order.

Later phases add Ng–Perron tests, response-surface p-values, bootstrap inference
(including bootstrap calibration of whole strategies) and outlier-robust tests.

## Installation

From GitHub:

```r
# install.packages("pak")
pak::pak("maarranz/UnitRoot-Tools-R")
```

## Specification

The econometric conventions (statistics, lag selection, critical values, the
result contract and the strategies) are defined in a language-neutral
specification, [`spec/specification.qmd`](spec/specification.qmd), shared with the
Python package UnitRoot-Tools. The two packages are independent implementations of
that specification and are validated against shared numerical fixtures.

## Licence

MIT © Miguel A. Arranz
