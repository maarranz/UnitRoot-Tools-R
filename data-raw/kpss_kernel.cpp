// KPSS statistic under the null of stationarity with iid N(0,1) errors and
// known long-run variance (= sample variance), for a large T: approximates
// the limiting distributions of eta_mu (int V1^2) and eta_tau (int V2^2).
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
NumericMatrix kpss_limit_stats(int T, int reps) {
  NumericMatrix out(reps, 2);
  std::vector<double> e(T);
  double st = 0, stt = 0; for (int t = 1; t <= T; ++t) { st += t; stt += (double)t * t; }
  double tbar = st / T, sxx = stt - T * tbar * tbar;
  RNGScope scope;
  for (int r = 0; r < reps; ++r) {
    NumericVector z = rnorm(T);
    double mean = 0; for (int t = 0; t < T; ++t) mean += z[t]; mean /= T;
    double sxy = 0; for (int t = 0; t < T; ++t) sxy += (t + 1 - tbar) * (z[t] - mean);
    double slope = sxy / sxx;
    for (int k = 0; k < 2; ++k) {
      double S = 0, sumS2 = 0, ss = 0;
      for (int t = 0; t < T; ++t) {
        double res = z[t] - mean - (k == 1 ? slope * (t + 1 - tbar) : 0.0);
        S += res; sumS2 += S * S; ss += res * res;
      }
      out(r, k) = sumS2 / ((double)T * T * (ss / T));
    }
  }
  return out;
}
