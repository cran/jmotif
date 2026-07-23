#include <jmotif.h>
#include <math.h>
// #include <iomanip>
// #include <iostream>
//
//' Computes a Piecewise Aggregate Approximation (PAA) for a time series.
//'
//' @param ts a timeseries to compute the PAA for.
//' @param paa_num the desired PAA size.
//' @useDynLib jmotif
//' @export
//' @references Keogh, E., Chakrabarti, K., Pazzani, M., Mehrotra, S.,
//' Dimensionality reduction for fast similarity search in large time series databases.
//' Knowledge and information Systems, 3(3), 263-286. (2001)
//' @examples
//' x = c(-1, -2, -1, 0, 2, 1, 1, 0)
//' x_paa3 = paa(x, 3)
//' #
//' plot(x, type = "l", main = c("8-points time series and its PAA transform into three points",
//'                           "PAA shown schematically in blue"))
//' points(x, pch = 16, lwd = 5)
//' #
//' paa_bounds = c(1, 1+7/3, 1+7/3*2, 8)
//' abline(v = paa_bounds, lty = 3, lwd = 2, col = "cornflowerblue")
//' segments(paa_bounds[1:3], x_paa3, paa_bounds[2:4], x_paa3, col = "cornflowerblue", lwd = 2)
//' points(x = c(1, 1+7/3, 1+7/3*2) + (7/3)/2, y = x_paa3, pch = 15, lwd = 5, col = "cornflowerblue")
// [[Rcpp::export]]
NumericVector paa(NumericVector ts, int paa_num) {
  return wrap(_paa2(Rcpp::as< std::vector<double> >(ts), paa_num));
}

std::vector<double> _paa2(const std::vector<double>& ts, int paa_num) {

  int len = ts.size();

  if (paa_num < 1) {
    stop("'paa_num' must be at least 1");
  }

  if(len < paa_num){
    stop("'paa_num' size is invalid");
  }

  if (len == paa_num) {
    return ts;
  }

  std::vector<double> res(paa_num, 0.0);
  double points_per_segment = (double) len / (double) paa_num;

  std::vector<double> breaks(paa_num + 1, 0);
  for(int i = 0; i < paa_num + 1; i++){
    breaks[i] = i * points_per_segment;
  }
  breaks[paa_num] = (double) len;

  for(int i = 0; i < paa_num; i++){

    double seg_start = breaks[i];
    double seg_end = breaks[i+1];

    double frac_begin = ceil(seg_start) - seg_start;
    double frac_end = seg_end - floor(seg_end);

    int full_begin = (int)floor(seg_start);
    int full_end = (int)ceil(seg_end);
    if(full_end > len) {
      full_end = len;
    }

    double sum_of_elems = 0.0;
    for (int j = full_begin; j < full_end; j++) {
      double v = ts[j];
      if (j == full_begin && frac_begin > 0) {
        v *= frac_begin;
      }
      if (j == full_end - 1 && frac_end > 0) {
        v *= frac_end;
      }
      sum_of_elems += v;
    }

    res[i] = sum_of_elems / points_per_segment;
  }
  return res;
}

void _paa2_range(const std::vector<double>& ts, int start, int end, int paa_num,
                 std::vector<double>& out) {

  int len = end - start;

  if(len < paa_num){
    stop("'paa_num' size is invalid");
  }

  if (len == paa_num) {
    out.resize(paa_num);
    for (int i = 0; i < paa_num; i++) {
      out[i] = ts[start + i];
    }
    return;
  }

  out.assign(paa_num, 0.0);
  double points_per_segment = (double) len / (double) paa_num;

  std::vector<double> breaks(paa_num + 1, 0);
  for(int i = 0; i < paa_num + 1; i++){
    breaks[i] = i * points_per_segment;
  }
  breaks[paa_num] = (double) len;

  for(int i = 0; i < paa_num; i++){

    double seg_start = breaks[i];
    double seg_end = breaks[i+1];

    double frac_begin = ceil(seg_start) - seg_start;
    double frac_end = seg_end - floor(seg_end);

    int full_begin = (int)floor(seg_start);
    int full_end = (int)ceil(seg_end);
    if(full_end > len) {
      full_end = len;
    }

    double sum_of_elems = 0.0;
    for (int j = full_begin; j < full_end; j++) {
      double v = ts[start + j];
      if (j == full_begin && frac_begin > 0) {
        v *= frac_begin;
      }
      if (j == full_end - 1 && frac_end > 0) {
        v *= frac_end;
      }
      sum_of_elems += v;
    }

    out[i] = sum_of_elems / points_per_segment;
  }
}
