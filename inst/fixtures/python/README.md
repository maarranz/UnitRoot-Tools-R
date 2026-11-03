# Reference files shared with the Python implementation

Byte-identical copies of two files of the Python package
[UnitRoot-Tools](https://github.com/maarranz/UnitRoot-Tools), taken from
commit 5a5a1ca (31 October 2026). Do not edit them: they are compared by
checksum in `tests/testthat/test-python_references.R`.

| File | Origin | SHA-256 |
|---|---|---|
| `adf_inputs.json` | `validation/scaled_svd/adf_inputs.json` | `0a5ef87c88547da466d37b1368e195037fe010de7a32a75e9234d8f2524ae8d4` |
| `scaled_references.json` | `validation/enders/scaled_references.json` | `a9f5ec22e573e42763f76fb43de03cbd75022a3e9bd5c99898f4d1e29ec1ce79` |

`adf_inputs.json` holds three series as exact binary64 values written in
hexadecimal: `y` (120 observations) and two copies rescaled to order
1e-12 (`positive`, and `negative` with the sign reversed).

`scaled_references.json` holds, for each series, trend `c` and `ct` and
0, 1 and 4 lagged differences (18 cases), the unrestricted ADF regression,
the regression restricted for phi1 (`c`) or phi3 (`ct`), and the auxiliary
regression with the unit root imposed. Coefficients, covariances, residual
sums of squares and t ratios were computed by exact binary64-to-decimal
elimination at 80 and 120 digits, with no production code involved
(`validation/enders/verify_scaled_references.py` in the Python repository
reproduces them). The file records the SHA-256 of the inputs it was built
from.

In these references the trend restarts at 1 in the regression sample,
whereas this package uses the original time index; only the constant (and
its t ratio) of the `ct` regressions depends on that, so it is not compared.
