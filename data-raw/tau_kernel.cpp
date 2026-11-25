// Null distribution of ADF-type t statistics, by sample size T and lag
// order k, for three deterministic cases. Written to be reused for the
// Engle-Granger test (cointegrationtools) through `nvars`:
//
//   nvars = 1 (unit-root test, unitroottools): y is a driftless Gaussian
//     random walk; the ADF regression of dy_t on [deterministics, y_{t-1},
//     dy_{t-1..t-k}] over the maximal sample t = k+2..T gives
//       case 0 "n":  no deterministics   (tau1)
//       case 1 "c":  constant            (tau2)
//       case 2 "ct": constant and trend  (tau3)
//
//   nvars = n > 1 (Engle-Granger): n independent driftless random walks
//     y_1..y_n; for each case, the cointegrating regression of y_1 on
//     [deterministics, y_2..y_n] over t = 1..T, and the ADF regression of its
//     residuals u_t on [u_{t-1}, du_{t-1..t-k}] over t = k+2..T, without
//     deterministics (they are in the cointegrating regression).
//
// Interface (R): tau_null_stats(T, lags, reps, nvars = 1) returns an array
// reps x length(lags) x 3, the third dimension being the cases n, c, ct.
// All lag orders use the same simulated paths. The RNG is R's (rnorm), so
// results are reproducible from R's seed (use L'Ecuyer streams in parallel).
#include <Rcpp.h>
#include <cmath>
#include <vector>
using namespace Rcpp;

struct TFit { double t; bool ok; };

// t ratio of the column at position `pos` of `cols` in the OLS regression on
// those columns, from the Gram system (G = X'X, g = X'y, yy = y'y; P columns
// in total, n observations). Cholesky on the subsystem.
static TFit t_ratio(const std::vector<double>& G, const std::vector<double>& g,
                    double yy, int P, int n, const std::vector<int>& cols, int pos) {
  int p = cols.size();
  TFit out{NA_REAL, false};
  if (n - p < 1) return out;
  std::vector<double> L(p * p, 0.0), z(p), beta(p), w(p, 0.0);
  for (int i = 0; i < p; ++i)
    for (int j = 0; j <= i; ++j) {
      double s = G[cols[i] * P + cols[j]];
      for (int m = 0; m < j; ++m) s -= L[i * p + m] * L[j * p + m];
      if (i == j) {
        if (!(s > 0)) return out;
        L[i * p + i] = std::sqrt(s);
      } else {
        L[i * p + j] = s / L[j * p + j];
      }
    }
  for (int i = 0; i < p; ++i) {
    double s = g[cols[i]];
    for (int m = 0; m < i; ++m) s -= L[i * p + m] * z[m];
    z[i] = s / L[i * p + i];
  }
  double fit_ss = 0.0;
  for (int i = 0; i < p; ++i) fit_ss += z[i] * z[i];
  double rss = yy - fit_ss;
  for (int i = p - 1; i >= 0; --i) {
    double s = z[i];
    for (int m = i + 1; m < p; ++m) s -= L[m * p + i] * beta[m];
    beta[i] = s / L[i * p + i];
  }
  for (int i = 0; i < p; ++i) {
    double s = (i == pos) ? 1.0 : 0.0;
    for (int m = 0; m < i; ++m) s -= L[i * p + m] * w[m];
    w[i] = s / L[i * p + i];
  }
  double v = 0.0;
  for (int i = 0; i < p; ++i) v += w[i] * w[i];
  double s2 = rss / (n - p);
  out.t = beta[pos] / std::sqrt(s2 * v);
  out.ok = true;
  return out;
}

// ADF t statistics of the series y (length T) for lag order k, for the
// requested column sets. Columns of the full design: 0 const, 1 trend (t of
// the response, 1-based), 2 y_{t-1}, 3.. lagged differences.
static void adf_taus(const std::vector<double>& y, int k, bool all_cases,
                     double* tau_n, double* tau_c, double* tau_ct) {
  int T = y.size(), n = T - 1 - k, P = 3 + k;
  std::vector<double> dy(T - 1);
  for (int t = 0; t < T - 1; ++t) dy[t] = y[t + 1] - y[t];
  std::vector<double> X(n * P), resp(n);
  for (int i = 0; i < n; ++i) {
    int idx = i + k;
    resp[i] = dy[idx];
    X[0 * n + i] = 1.0;
    X[1 * n + i] = idx + 2;
    X[2 * n + i] = y[idx];
    for (int l = 1; l <= k; ++l) X[(2 + l) * n + i] = dy[idx - l];
  }
  std::vector<double> G(P * P), g(P);
  double yy = 0.0;
  for (int i = 0; i < n; ++i) yy += resp[i] * resp[i];
  for (int a = 0; a < P; ++a) {
    for (int b = 0; b <= a; ++b) {
      double s = 0;
      for (int i = 0; i < n; ++i) s += X[a * n + i] * X[b * n + i];
      G[a * P + b] = s; G[b * P + a] = s;
    }
    double s = 0;
    for (int i = 0; i < n; ++i) s += X[a * n + i] * resp[i];
    g[a] = s;
  }
  std::vector<int> lagcols;
  for (int l = 0; l < k; ++l) lagcols.push_back(3 + l);
  std::vector<int> cn = {2};       cn.insert(cn.end(), lagcols.begin(), lagcols.end());
  *tau_n = t_ratio(G, g, yy, P, n, cn, 0).t;
  if (all_cases) {
    std::vector<int> cc = {0, 2};    cc.insert(cc.end(), lagcols.begin(), lagcols.end());
    std::vector<int> cct = {0, 1, 2}; cct.insert(cct.end(), lagcols.begin(), lagcols.end());
    *tau_c = t_ratio(G, g, yy, P, n, cc, 1).t;
    *tau_ct = t_ratio(G, g, yy, P, n, cct, 2).t;
  }
}

// Residuals of the regression of y1 on [deterministics (case), y2..yn] over
// t = 1..T (normal equations; small number of regressors).
static std::vector<double> coint_residuals(const std::vector<std::vector<double>>& Y,
                                           int cs) {
  int T = Y[0].size(), nv = Y.size();
  int d = (cs == 0) ? 0 : (cs == 1 ? 1 : 2);
  int p = d + nv - 1;
  std::vector<double> Xc(T * p);
  for (int t = 0; t < T; ++t) {
    int c = 0;
    if (d >= 1) Xc[(c++) * T + t] = 1.0;
    if (d >= 2) Xc[(c++) * T + t] = t + 1;
    for (int v = 1; v < nv; ++v) Xc[(c++) * T + t] = Y[v][t];
  }
  std::vector<double> G(p * p), g(p);
  for (int a = 0; a < p; ++a) {
    for (int b = 0; b <= a; ++b) {
      double s = 0;
      for (int t = 0; t < T; ++t) s += Xc[a * T + t] * Xc[b * T + t];
      G[a * p + b] = s; G[b * p + a] = s;
    }
    double s = 0;
    for (int t = 0; t < T; ++t) s += Xc[a * T + t] * Y[0][t];
    g[a] = s;
  }
  // Cholesky solve
  std::vector<double> L(p * p, 0.0), z(p), beta(p);
  for (int i = 0; i < p; ++i)
    for (int j = 0; j <= i; ++j) {
      double s = G[i * p + j];
      for (int m = 0; m < j; ++m) s -= L[i * p + m] * L[j * p + m];
      if (i == j) L[i * p + i] = std::sqrt(s); else L[i * p + j] = s / L[j * p + j];
    }
  for (int i = 0; i < p; ++i) {
    double s = g[i];
    for (int m = 0; m < i; ++m) s -= L[i * p + m] * z[m];
    z[i] = s / L[i * p + i];
  }
  for (int i = p - 1; i >= 0; --i) {
    double s = z[i];
    for (int m = i + 1; m < p; ++m) s -= L[m * p + i] * beta[m];
    beta[i] = s / L[i * p + i];
  }
  std::vector<double> u(T);
  for (int t = 0; t < T; ++t) {
    double f = 0;
    for (int a = 0; a < p; ++a) f += Xc[a * T + t] * beta[a];
    u[t] = Y[0][t] - f;
  }
  return u;
}

// [[Rcpp::export]]
NumericVector tau_null_stats(int T, IntegerVector lags, int reps, int nvars = 1) {
  if (nvars < 1) stop("nvars must be at least 1");
  int nk = lags.size();
  NumericVector out(reps * nk * 3);
  out.attr("dim") = IntegerVector::create(reps, nk, 3);
  RNGScope scope;
  std::vector<std::vector<double>> Y(nvars, std::vector<double>(T));
  for (int r = 0; r < reps; ++r) {
    for (int v = 0; v < nvars; ++v) {
      NumericVector e = rnorm(T);
      double acc = 0.0;
      for (int t = 0; t < T; ++t) { acc += e[t]; Y[v][t] = acc; }
    }
    if (nvars == 1) {
      for (int j = 0; j < nk; ++j) {
        double a, b, c;
        adf_taus(Y[0], lags[j], true, &a, &b, &c);
        out[r + reps * (j + nk * 0)] = a;
        out[r + reps * (j + nk * 1)] = b;
        out[r + reps * (j + nk * 2)] = c;
      }
    } else {
      for (int cs = 0; cs < 3; ++cs) {
        std::vector<double> u = coint_residuals(Y, cs);
        for (int j = 0; j < nk; ++j) {
          double a, unused1, unused2;
          adf_taus(u, lags[j], false, &a, &unused1, &unused2);
          out[r + reps * (j + nk * cs)] = a;
        }
      }
    }
  }
  return out;
}
