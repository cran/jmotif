#include <jmotif.h>

//' Z-normalizes a time series by subtracting its mean and dividing by the standard deviation.
//'
//' @param ts the input time series.
//' @param threshold the z-normalization threshold value, if the input time series' standard
//' deviation will be found less than this value, the procedure will not be applied,
//' so the "under-threshold-noise" would not get amplified.
//' @useDynLib jmotif
//' @export
//' @references Dina Goldin and Paris Kanellakis,
//' On similarity queries for time-series data: Constraint specification and implementation.
//' In Principles and Practice of Constraint Programming (CP 1995), pages 137-153. (1995)
//' @examples
//' x = seq(0, pi*4, 0.02)
//' y = sin(x) * 5 + rnorm(length(x))
//' plot(x, y, type="l", col="blue")
//' lines(x, znorm(y, 0.01), type="l", col="red")
// [[Rcpp::export]]
NumericVector znorm(NumericVector ts, double threshold = 0.01) {

  if (ts.size() == 0) stop("Input vector cannot be empty");

  std::vector<double> ts_std = as<std::vector<double>>(ts);
  std::vector<double> out = _znorm(ts_std, threshold);
  return wrap(out);
}

std::vector<double> _znorm(const std::vector<double>& ts, double threshold) {

  double sum = std::accumulate(ts.begin(), ts.end(), 0.0);
  double mean =  sum / ts.size();
  // Rcout << " mean2 " << mean << "\n";

  std::vector<double> diff(ts.size());
  //std::transform(ts.begin(), ts.end(), diff.begin(), std::bind2nd(std::minus<double>(), mean));
  for(unsigned i=0; i<ts.size(); i++)
    diff[i] = ts[i]-mean;

  double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
  // Population standard deviation (divide by n, not n-1). This matches the
  // Matrix Profile / MASS convention and the saxpy reference, and makes each
  // normalized window have empirical variance exactly 1 -- the assumption
  // behind SAX's equiprobable Gaussian breakpoints. (Aligned 2026-06; was n-1.)
  double stdev = std::sqrt(sq_sum / ts.size());
  // Rcout << " stdev2 " << stdev << "\n";

  if (stdev < threshold) return ts; // return clone

  for (unsigned i = 0; i < ts.size(); i++)
    diff[i] = (ts[i] - mean) / stdev;

  return diff;

}

void _znorm_slice(const std::vector<double>& ts, int start, int end,
                  double threshold, std::vector<double>& out) {
  int len = end - start;
  out.resize(len);

  double sum = 0.0;
  for (int i = start; i < end; i++) {
    sum += ts[i];
  }
  double mean = sum / len;

  double sq_sum = 0.0;
  for (int i = start; i < end; i++) {
    double d = ts[i] - mean;
    sq_sum += d * d;
  }
  double stdev = std::sqrt(sq_sum / len);

  if (stdev < threshold) {
    for (int i = start; i < end; i++) {
      out[i - start] = ts[i];
    }
    return;
  }

  for (int i = start; i < end; i++) {
    out[i - start] = (ts[i] - mean) / stdev;
  }
}
