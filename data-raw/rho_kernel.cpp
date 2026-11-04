// Dickey-Fuller normalized-bias statistic n * b_hat (n = T - 1) under the
// unit-root null (random walk, N(0,1)), for the regressions with a constant
// and with a constant and trend; no lagged differences. Its limit is that of
// Phillips-Perron's Z_alpha.
#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;
// [[Rcpp::export]]
NumericMatrix rho_null_stats(int T, int reps) {
  NumericMatrix out(reps, 2);
  int n = T - 1;
  RNGScope scope;
  std::vector<double> y(T);
  for (int r = 0; r < reps; ++r) {
    NumericVector e = rnorm(T);
    double acc = 0; for (int t = 0; t < T; ++t) { acc += e[t]; y[t] = acc; }
    // regressors: 1, t (= 2..T), y_{t-1}; response dy_t
    double S1 = n, St = 0, Stt = 0, Sy = 0, Syy = 0, Sty = 0, Sd = 0, Sdy = 0, Sdt = 0;
    for (int i = 0; i < n; ++i) {
      double t = i + 2, yl = y[i], d = y[i + 1] - y[i];
      St += t; Stt += t * t; Sy += yl; Syy += yl * yl; Sty += t * yl; Sd += d; Sdy += d * yl; Sdt += d * t;
    }
    // constant: b = cov(d, yl) / var(yl)
    double vyy = Syy - Sy * Sy / S1, cdy = Sdy - Sd * Sy / S1;
    out(r, 0) = n * cdy / vyy;
    // constant and trend: partial out (1, t) by Frisch-Waugh
    double vtt = Stt - St * St / S1, cty = Sty - St * Sy / S1, cdt = Sdt - Sd * St / S1;
    double vyy_t = vyy - cty * cty / vtt, cdy_t = cdy - cdt * cty / vtt;
    out(r, 1) = n * cdy_t / vyy_t;
  }
  return out;
}
