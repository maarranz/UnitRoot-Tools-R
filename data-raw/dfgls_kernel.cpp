// DF-GLS statistics under the null (Elliott, Rothenberg and Stock 1996), for
// all requested lag orders on the same simulated path.
#include <Rcpp.h>
#include <cmath>
#include <vector>
using namespace Rcpp;

// OLS t ratio of the first column of X (n x p, column major) on y, by
// Cholesky of X'X (p <= 13, well conditioned for these designs).
static double t_first(const std::vector<double>& X, const std::vector<double>& y,
                      int n, int p) {
  std::vector<double> A(p * p, 0.0), b(p, 0.0);
  for (int i = 0; i < p; ++i) {
    for (int j = 0; j <= i; ++j) {
      double s = 0.0;
      const double* xi = &X[i * n]; const double* xj = &X[j * n];
      for (int t = 0; t < n; ++t) s += xi[t] * xj[t];
      A[i * p + j] = s; A[j * p + i] = s;
    }
    double s = 0.0; const double* xi = &X[i * n];
    for (int t = 0; t < n; ++t) s += xi[t] * y[t];
    b[i] = s;
  }
  // Cholesky A = L L'
  std::vector<double> L(p * p, 0.0);
  for (int i = 0; i < p; ++i) {
    for (int j = 0; j <= i; ++j) {
      double s = A[i * p + j];
      for (int k = 0; k < j; ++k) s -= L[i * p + k] * L[j * p + k];
      if (i == j) L[i * p + i] = std::sqrt(s); else L[i * p + j] = s / L[j * p + j];
    }
  }
  // solve for beta
  std::vector<double> z(p), beta(p);
  for (int i = 0; i < p; ++i) { double s = b[i]; for (int k = 0; k < i; ++k) s -= L[i * p + k] * z[k]; z[i] = s / L[i * p + i]; }
  for (int i = p - 1; i >= 0; --i) { double s = z[i]; for (int k = i + 1; k < p; ++k) s -= L[k * p + i] * beta[k]; beta[i] = s / L[i * p + i]; }
  // (X'X)^{-1}_{11}: solve L w = e1, then ||w||^2 over the inverse of L' ... use
  // e1' A^{-1} e1 = || L^{-1} e1 ||^2
  std::vector<double> w(p, 0.0);
  for (int i = 0; i < p; ++i) { double s = (i == 0) ? 1.0 : 0.0; for (int k = 0; k < i; ++k) s -= L[i * p + k] * w[k]; w[i] = s / L[i * p + i]; }
  double inv11 = 0.0; for (int i = 0; i < p; ++i) inv11 += w[i] * w[i];
  double ssr = 0.0;
  for (int t = 0; t < n; ++t) {
    double f = 0.0; for (int i = 0; i < p; ++i) f += X[i * n + t] * beta[i];
    double e = y[t] - f; ssr += e * e;
  }
  double s2 = ssr / (n - p);
  return beta[0] / std::sqrt(s2 * inv11);
}

// [[Rcpp::export]]
NumericMatrix dfgls_null_stats(int T, IntegerVector lags, bool trend, int reps) {
  int nk = lags.size();
  NumericMatrix out(reps, nk);
  double cbar = trend ? -13.5 : -7.0, a = 1.0 + cbar / T;
  std::vector<double> y(T), yd(T), dy(T - 1);
  RNGScope scope;
  for (int r = 0; r < reps; ++r) {
    NumericVector e = rnorm(T);
    double acc = 0.0;
    for (int t = 0; t < T; ++t) { acc += e[t]; y[t] = acc; }
    // GLS detrending: regress quasi-differenced y on quasi-differenced z
    if (!trend) {
      double num = y[0], den = 1.0, d = 1.0 - a;
      for (int t = 1; t < T; ++t) { num += d * (y[t] - a * y[t - 1]); den += d * d; }
      double mu = num / den;
      for (int t = 0; t < T; ++t) yd[t] = y[t] - mu;
    } else {
      // z_t = (1, t), t = 1..T
      double s11 = 1.0, s12 = 1.0, s22 = 1.0, r1 = y[0], r2 = y[0];
      for (int t = 1; t < T; ++t) {
        double z1 = 1.0 - a, z2 = (t + 1) - a * t, q = y[t] - a * y[t - 1];
        s11 += z1 * z1; s12 += z1 * z2; s22 += z2 * z2; r1 += z1 * q; r2 += z2 * q;
      }
      double det = s11 * s22 - s12 * s12;
      double b1 = (s22 * r1 - s12 * r2) / det, b2 = (s11 * r2 - s12 * r1) / det;
      for (int t = 0; t < T; ++t) yd[t] = y[t] - b1 - b2 * (t + 1);
    }
    for (int t = 0; t < T - 1; ++t) dy[t] = yd[t + 1] - yd[t];
    for (int j = 0; j < nk; ++j) {
      int k = lags[j], n = T - 1 - k, p = 1 + k;
      std::vector<double> X(n * p), resp(n);
      for (int t = 0; t < n; ++t) {
        int idx = t + k;            // dy index of the response
        resp[t] = dy[idx];
        X[t] = yd[idx];             // yd_{t-1}
        for (int l = 1; l <= k; ++l) X[l * n + t] = dy[idx - l];
      }
      out(r, j) = t_first(X, resp, n, p);
    }
  }
  return out;
}
