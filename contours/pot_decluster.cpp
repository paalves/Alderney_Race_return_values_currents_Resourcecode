#include <Rcpp.h>
#include <deque>
#include <vector>
#include <algorithm>
#include <cmath>
#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp17)]]

using namespace Rcpp;

// Type-7 sample quantile of a sorted vector.
static inline double q7(const std::vector<double>& x, double p) {
  const std::size_t n = x.size();
  if (!n) return NA_REAL;
  if (n == 1 || p <= 0) return x[0];
  if (p >= 1) return x[n - 1];
  const double h = (n - 1) * p;
  const std::size_t lo = std::floor(h), hi = lo + 1;
  if (hi >= n) return x[n - 1];
  return x[lo] + (h - lo) * (x[hi] - x[lo]);
}

struct DirResult {
  std::vector<double> peaks, exceedances;
  double mu = NA_REAL;
  int n_peaks = 0;
};

// Project the data on one direction, keep moving-window local maxima as peaks and return the excesses above their (1 - zeta) quantile.
static DirResult process_direction(const double* Y0, const double* Y1, const double* Y2, int n,
                                   double d0, double d1, double d2, int h, double zeta) {
  DirResult out;
  const int k = 2 * h + 1, half = (k - 1) / 2;
  if (k > n || k < 1 || n <= 0) return out;

  std::vector<double> r(n), rmax(n, NA_REAL);
  std::deque<int> dq;
  for (int i = 0; i < n; ++i) r[i] = Y0[i] * d0 + Y1[i] * d1 + Y2[i] * d2;
  for (int i = 0; i < k - 1; ++i) {
    while (!dq.empty() && r[dq.back()] <= r[i]) dq.pop_back();
    dq.push_back(i);
  }
  for (int i = k - 1; i < n; ++i) {
    while (!dq.empty() && r[dq.back()] <= r[i]) dq.pop_back();
    dq.push_back(i);
    while (!dq.empty() && dq.front() <= i - k) dq.pop_front();
    rmax[i - half] = r[dq.front()];
  }

  std::vector<int> idx;
  idx.reserve(n / k + 16);
  for (int i = 0; i < n; ++i)
    if (!std::isnan(rmax[i]) && std::fabs(r[i] - rmax[i]) < 1e-12) idx.push_back(i);

  if (idx.size() > 1) {
    std::vector<int> keep{idx[0]};
    for (std::size_t i = 1; i < idx.size(); ++i)
      if (idx[i] - keep.back() >= h) keep.push_back(idx[i]);
    idx.swap(keep);
  }

  out.n_peaks = idx.size();
  out.peaks.reserve(idx.size());
  for (int i : idx) out.peaks.push_back(r[i]);
  if (out.peaks.empty()) return out;

  std::vector<double> s = out.peaks;
  std::sort(s.begin(), s.end());
  out.mu = q7(s, 1 - zeta);
  for (double x : out.peaks) if (x > out.mu) out.exceedances.push_back(x - out.mu);
  return out;
}

// Decluster the projections on every direction in parallel and return peaks, thresholds and excesses.
// [[Rcpp::export]]
List pot_decluster_all(NumericMatrix Y, NumericMatrix dirs, int h_hours, double zeta, int n_threads = 0) {
  if (Y.ncol() != 3) stop("Y must have 3 columns");
  if (dirs.ncol() != 3) stop("dirs must have 3 columns");

  const int n = Y.nrow(), nd = dirs.nrow();
  const double *Y0 = &Y(0, 0), *Y1 = &Y(0, 1), *Y2 = &Y(0, 2);
  std::vector<double> d(3 * nd);
  for (int j = 0; j < nd; ++j) {
    d[3 * j] = dirs(j, 0);
    d[3 * j + 1] = dirs(j, 1);
    d[3 * j + 2] = dirs(j, 2);
  }

  std::vector<DirResult> res(nd);
#ifdef _OPENMP
  if (n_threads > 0) omp_set_num_threads(n_threads);
  #pragma omp parallel for schedule(dynamic)
#endif
  for (int j = 0; j < nd; ++j)
    res[j] = process_direction(Y0, Y1, Y2, n, d[3 * j], d[3 * j + 1], d[3 * j + 2], h_hours, zeta);

  List peaks(nd), exc(nd);
  NumericVector mu(nd);
  IntegerVector n_peaks(nd);
  for (int j = 0; j < nd; ++j) {
    peaks[j] = NumericVector(res[j].peaks.begin(), res[j].peaks.end());
    exc[j] = NumericVector(res[j].exceedances.begin(), res[j].exceedances.end());
    mu[j] = res[j].mu;
    n_peaks[j] = res[j].n_peaks;
  }
  return List::create(_["peaks"] = peaks, _["exceedances"] = exc, _["mu"] = mu, _["n_peaks"] = n_peaks);
}
