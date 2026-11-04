// phi1 (constant), phi2 and phi3 (constant and trend) statistics of the ADF
// regression under the unit-root null, with tau2 and tau3 as a check, for
// all requested lag orders on the same simulated random walk. All five
// regressions are column subsets of [const, trend, y_{t-1}, dy_{t-1..t-k}]
// on the maximal sample t = k+2..T, so one cross-product matrix serves all.
#include <Rcpp.h>
#include <cmath>
#include <vector>
using namespace Rcpp;

struct Fit { double rss; double b0; double v00; };

// OLS on the columns `cols` of the Gram system (G = X'X, g = X'y, yy = y'y).
// Returns RSS and, for the column at position `pos` in `cols`, the
// coefficient and the (X'X)^{-1} diagonal element.
static Fit sub_fit(const std::vector<double>& G, const std::vector<double>& g,
                   double yy, int P, const std::vector<int>& cols, int pos) {
  int p = cols.size();
  Fit f{yy, 0.0, 0.0};
  if (p == 0) return f;
  std::vector<double> L(p * p, 0.0), z(p), beta(p);
  for (int i = 0; i < p; ++i)
    for (int j = 0; j <= i; ++j) {
      double s = G[cols[i] * P + cols[j]];
      for (int k = 0; k < j; ++k) s -= L[i * p + k] * L[j * p + k];
      if (i == j) L[i * p + i] = std::sqrt(s); else L[i * p + j] = s / L[j * p + j];
    }
  for (int i = 0; i < p; ++i) { double s = g[cols[i]]; for (int k = 0; k < i; ++k) s -= L[i * p + k] * z[k]; z[i] = s / L[i * p + i]; }
  double fit_ss = 0.0; for (int i = 0; i < p; ++i) fit_ss += z[i] * z[i];
  f.rss = yy - fit_ss;
  if (pos >= 0) {
    for (int i = p - 1; i >= 0; --i) { double s = z[i]; for (int k = i + 1; k < p; ++k) s -= L[k * p + i] * beta[k]; beta[i] = s / L[i * p + i]; }
    std::vector<double> w(p, 0.0);
    for (int i = 0; i < p; ++i) { double s = (i == pos) ? 1.0 : 0.0; for (int k = 0; k < i; ++k) s -= L[i * p + k] * w[k]; w[i] = s / L[i * p + i]; }
    // e_pos' A^{-1} e_pos = || L^{-1} e_pos ||^2 (w is zero above pos)
    double v = 0.0; for (int i = 0; i < p; ++i) v += w[i] * w[i];
    f.b0 = beta[pos]; f.v00 = v;
  }
  return f;
}

// [[Rcpp::export]]
NumericVector phi_null_stats(int T, IntegerVector lags, int reps) {
  int nk = lags.size();
  // dims: reps x nk x 5 (phi1, phi2, phi3, tau2, tau3)
  NumericVector out(reps * nk * 5);
  out.attr("dim") = IntegerVector::create(reps, nk, 5);
  std::vector<double> y(T), dy(T - 1);
  RNGScope scope;
  for (int r = 0; r < reps; ++r) {
    NumericVector e = rnorm(T);
    double acc = 0.0;
    for (int t = 0; t < T; ++t) { acc += e[t]; y[t] = acc; }
    for (int t = 0; t < T - 1; ++t) dy[t] = y[t + 1] - y[t];
    for (int j = 0; j < nk; ++j) {
      int k = lags[j], n = T - 1 - k, P = 3 + k;
      // columns: 0 const, 1 trend (original index t = k+2..T), 2 y_{t-1}, 3.. lags
      std::vector<double> X(n * P), resp(n);
      for (int i = 0; i < n; ++i) {
        int idx = i + k;                         // dy index of the response
        resp[i] = dy[idx];
        X[0 * n + i] = 1.0;
        X[1 * n + i] = idx + 2;                  // t (1-based) of the response
        X[2 * n + i] = y[idx];                   // y_{t-1}
        for (int l = 1; l <= k; ++l) X[(2 + l) * n + i] = dy[idx - l];
      }
      std::vector<double> G(P * P), g(P);
      double yy = 0.0; for (int i = 0; i < n; ++i) yy += resp[i] * resp[i];
      for (int a = 0; a < P; ++a) {
        for (int b = 0; b <= a; ++b) { double s = 0; for (int i = 0; i < n; ++i) s += X[a * n + i] * X[b * n + i]; G[a * P + b] = s; G[b * P + a] = s; }
        double s = 0; for (int i = 0; i < n; ++i) s += X[a * n + i] * resp[i]; g[a] = s;
      }
      std::vector<int> lagcols; for (int l = 0; l < k; ++l) lagcols.push_back(3 + l);
      std::vector<int> c_u = {0, 2}; c_u.insert(c_u.end(), lagcols.begin(), lagcols.end());
      std::vector<int> ct_u = {0, 1, 2}; ct_u.insert(ct_u.end(), lagcols.begin(), lagcols.end());
      std::vector<int> c0 = {0}; c0.insert(c0.end(), lagcols.begin(), lagcols.end());
      Fit fc = sub_fit(G, g, yy, P, c_u, 1);
      Fit fct = sub_fit(G, g, yy, P, ct_u, 2);
      Fit fl = sub_fit(G, g, yy, P, lagcols, -1);
      Fit fc0 = sub_fit(G, g, yy, P, c0, -1);
      double s2c = fc.rss / (n - (int)c_u.size()), s2ct = fct.rss / (n - (int)ct_u.size());
      double phi1 = ((fl.rss - fc.rss) / 2.0) / s2c;
      double phi2 = ((fl.rss - fct.rss) / 3.0) / s2ct;
      double phi3 = ((fc0.rss - fct.rss) / 2.0) / s2ct;
      double tau2 = fc.b0 / std::sqrt(s2c * fc.v00), tau3 = fct.b0 / std::sqrt(s2ct * fct.v00);
      double vals[5] = {phi1, phi2, phi3, tau2, tau3};
      for (int s = 0; s < 5; ++s) out[r + reps * (j + nk * s)] = vals[s];
    }
  }
  return out;
}
