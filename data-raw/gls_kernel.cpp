// GLS-based unit-root statistics under the null (random walk, N(0,1)):
// DF-GLS t, Ng-Perron MZa, MZt, MSB, MPT and the ERS point-optimal PT, for
// all requested lag orders on the same path. Definitions as in
// R/gls_tests.R (gls_statistics) and dfgls().
#include <Rcpp.h>
#include <cmath>
#include <vector>
using namespace Rcpp;

// OLS of y on X (n x p, column major) by Cholesky of X'X; returns SSR, the
// coefficients and (X'X)^{-1}_{11}.
static void ols(const std::vector<double>& X, const std::vector<double>& y, int n, int p,
                double& ssr, std::vector<double>& beta, double& inv11) {
  std::vector<double> A(p * p), b(p), L(p * p, 0.0), z(p);
  beta.assign(p, 0.0);
  for (int i = 0; i < p; ++i) {
    for (int j = 0; j <= i; ++j) { double s = 0; for (int t = 0; t < n; ++t) s += X[i*n+t] * X[j*n+t]; A[i*p+j] = s; A[j*p+i] = s; }
    double s = 0; for (int t = 0; t < n; ++t) s += X[i*n+t] * y[t]; b[i] = s;
  }
  for (int i = 0; i < p; ++i) for (int j = 0; j <= i; ++j) {
    double s = A[i*p+j]; for (int k = 0; k < j; ++k) s -= L[i*p+k] * L[j*p+k];
    if (i == j) L[i*p+i] = std::sqrt(s); else L[i*p+j] = s / L[j*p+j];
  }
  for (int i = 0; i < p; ++i) { double s = b[i]; for (int k = 0; k < i; ++k) s -= L[i*p+k] * z[k]; z[i] = s / L[i*p+i]; }
  for (int i = p - 1; i >= 0; --i) { double s = z[i]; for (int k = i + 1; k < p; ++k) s -= L[k*p+i] * beta[k]; beta[i] = s / L[i*p+i]; }
  std::vector<double> w(p, 0.0);
  for (int i = 0; i < p; ++i) { double s = (i == 0) ? 1.0 : 0.0; for (int k = 0; k < i; ++k) s -= L[i*p+k] * w[k]; w[i] = s / L[i*p+i]; }
  inv11 = 0; for (int i = 0; i < p; ++i) inv11 += w[i] * w[i];
  ssr = 0;
  for (int t = 0; t < n; ++t) { double f = 0; for (int i = 0; i < p; ++i) f += X[i*n+t] * beta[i]; double e = y[t] - f; ssr += e * e; }
}

// GLS regression of the quasi-differenced y on the quasi-differenced z
// (z = 1 or (1, t)); returns SSR and the coefficients.
static void gls_fit(const std::vector<double>& y, int T, bool trend, double a,
                    double& ssr, double& b1, double& b2) {
  int p = trend ? 2 : 1;
  std::vector<double> X(T * p), q(T), beta; double inv11;
  q[0] = y[0]; X[0] = 1.0; if (trend) X[T] = 1.0;
  for (int t = 1; t < T; ++t) {
    q[t] = y[t] - a * y[t - 1];
    X[t] = 1.0 - a;
    if (trend) X[T + t] = (t + 1) - a * t;
  }
  ols(X, q, T, p, ssr, beta, inv11);
  b1 = beta[0]; b2 = trend ? beta[1] : 0.0;
}

// [[Rcpp::export]]
NumericVector gls_null_stats(int T, IntegerVector lags, bool trend, int reps) {
  int nk = lags.size(), ns = 6;   // DFGLS, MZa, MZt, MSB, MPT, PT
  NumericVector out(reps * nk * ns);
  out.attr("dim") = IntegerVector::create(reps, nk, ns);
  double cbar = trend ? -13.5 : -7.0, abar = 1.0 + cbar / T;
  std::vector<double> y(T), yd(T), dy(T - 1);
  RNGScope scope;
  for (int r = 0; r < reps; ++r) {
    NumericVector e = rnorm(T);
    double acc = 0; for (int t = 0; t < T; ++t) { acc += e[t]; y[t] = acc; }
    double s_abar, s_one, b1, b2, d1, d2;
    gls_fit(y, T, trend, abar, s_abar, b1, b2);
    gls_fit(y, T, trend, 1.0, s_one, d1, d2);
    for (int t = 0; t < T; ++t) yd[t] = y[t] - b1 - b2 * (t + 1);
    for (int t = 0; t < T - 1; ++t) dy[t] = yd[t + 1] - yd[t];
    double kappa = 0; for (int t = 0; t < T - 1; ++t) kappa += yd[t] * yd[t]; kappa /= (double)T * T;
    double last = yd[T - 1] * yd[T - 1] / T;
    for (int j = 0; j < nk; ++j) {
      int k = lags[j], n = T - 1 - k, p = 1 + k;
      std::vector<double> X(n * p), resp(n), beta; double ssr, inv11;
      for (int i = 0; i < n; ++i) {
        int idx = i + k; resp[i] = dy[idx]; X[i] = yd[idx];
        for (int l = 1; l <= k; ++l) X[l * n + i] = dy[idx - l];
      }
      ols(X, resp, n, p, ssr, beta, inv11);
      double tdf = beta[0] / std::sqrt(ssr / (n - p) * inv11);
      double sb = 0; for (int l = 1; l <= k; ++l) sb += beta[l];
      double s2 = (ssr / n) / ((1 - sb) * (1 - sb));
      double mza = (last - s2) / (2 * kappa), msb = std::sqrt(kappa / s2);
      double mpt = trend ? (cbar * cbar * kappa + (1 - cbar) * last) / s2
                         : (cbar * cbar * kappa - cbar * last) / s2;
      double pt = (s_abar - abar * s_one) / s2;
      double v[6] = {tdf, mza, mza * msb, msb, mpt, pt};
      for (int s = 0; s < ns; ++s) out[r + reps * (j + nk * s)] = v[s];
    }
  }
  return out;
}
