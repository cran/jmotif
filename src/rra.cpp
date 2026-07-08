#include <jmotif.h>

#include <limits>
#include <chrono>
#include <ctime>
#include <vector>
#include <random>
#include <algorithm> // std::stable_sort, std::max


class rule_interval {
public:
  int rule_id;
  int start;
  int end;
  double cover;
};

struct sort_intervals {
  bool operator()(const rule_interval &left, const rule_interval &right) {
    return (left.cover < right.cover);
  }
};

// Per-point normalized Euclidean distance between two subsequences.
//
// Both subsequences are Z-NORMALIZED before the distance is taken, so the
// comparison is of shape, not amplitude/offset -- matching the HOT-SAX and
// brute-force discord routines (hot-sax.cpp / discord.cpp both _znorm their
// windows) and the canonical jMotif GrammarViz RRA. When the spans differ in
// length, the longer is first PAA-reduced (_paa2) to the shorter length. The
// result is euclidean(znorm(a), znorm(b)) / count, divided by the number of
// compared points so spans of different lengths stay comparable.
namespace {

double slice_mean(const std::vector<double>& ts, int start, int end) {
  double sum = 0.0;
  for (int i = start; i < end; i++) {
    sum += ts[i];
  }
  return sum / (end - start);
}

double slice_stdev_pop(const std::vector<double>& ts, int start, int end, double mean) {
  double sq_sum = 0.0;
  for (int i = start; i < end; i++) {
    double d = ts[i] - mean;
    sq_sum += d * d;
  }
  return std::sqrt(sq_sum / (end - start));
}

double vec_mean(const std::vector<double>& v) {
  double sum = 0.0;
  for (double x : v) {
    sum += x;
  }
  return sum / v.size();
}

double vec_stdev_pop(const std::vector<double>& v, double mean) {
  double sq_sum = 0.0;
  for (double x : v) {
    double d = x - mean;
    sq_sum += d * d;
  }
  return std::sqrt(sq_sum / v.size());
}

double znormed_rms_dist_equal(
    const std::vector<double>& ts, int s1, int e1, int s2, int e2,
    double n_threshold) {
  int n = e1 - s1;
  double m1 = slice_mean(ts, s1, e1);
  double m2 = slice_mean(ts, s2, e2);
  double sd1 = slice_stdev_pop(ts, s1, e1, m1);
  double sd2 = slice_stdev_pop(ts, s2, e2, m2);
  bool raw1 = sd1 < n_threshold;
  bool raw2 = sd2 < n_threshold;

  double res = 0.0;
  for (int i = 0; i < n; i++) {
    double a = raw1 ? ts[s1 + i] : (ts[s1 + i] - m1) / sd1;
    double b = raw2 ? ts[s2 + i] : (ts[s2 + i] - m2) / sd2;
    double d = a - b;
    res += d * d;
  }
  return std::sqrt(res) / n;
}

double znormed_rms_dist_slice_vec(
    const std::vector<double>& ts, int s1, int e1,
    const std::vector<double>& b, double n_threshold) {
  int n = e1 - s1;
  double m1 = slice_mean(ts, s1, e1);
  double m2 = vec_mean(b);
  double sd1 = slice_stdev_pop(ts, s1, e1, m1);
  double sd2 = vec_stdev_pop(b, m2);
  bool raw1 = sd1 < n_threshold;
  bool raw2 = sd2 < n_threshold;

  double res = 0.0;
  for (int i = 0; i < n; i++) {
    double a = raw1 ? ts[s1 + i] : (ts[s1 + i] - m1) / sd1;
    double bv = raw2 ? b[i] : (b[i] - m2) / sd2;
    double d = a - bv;
    res += d * d;
  }
  return std::sqrt(res) / n;
}

double znormed_rms_dist_vec_slice(
    const std::vector<double>& a,
    const std::vector<double>& ts, int s2, int e2,
    double n_threshold) {
  int n = a.size();
  double m1 = vec_mean(a);
  double m2 = slice_mean(ts, s2, e2);
  double sd1 = vec_stdev_pop(a, m1);
  double sd2 = slice_stdev_pop(ts, s2, e2, m2);
  bool raw1 = sd1 < n_threshold;
  bool raw2 = sd2 < n_threshold;

  double res = 0.0;
  for (int i = 0; i < n; i++) {
    double av = raw1 ? a[i] : (a[i] - m1) / sd1;
    double b = raw2 ? ts[s2 + i] : (ts[s2 + i] - m2) / sd2;
    double d = av - b;
    res += d * d;
  }
  return std::sqrt(res) / n;
}

double znormed_rms_dist_precomputed(
    const std::vector<double>& za, const std::vector<double>& zb) {
  double res = 0.0;
  int n = za.size();
  for (int i = 0; i < n; i++) {
    double d = za[i] - zb[i];
    res += d * d;
  }
  return std::sqrt(res) / n;
}

} // namespace

double _normalized_distance(int start1, int end1, int start2, int end2,
                            const std::vector<double>& series, double n_threshold,
                            int w_size,
                            const std::vector<std::vector<double>>* znorm_windows) {

  int len1 = end1 - start1;
  int len2 = end2 - start2;

  if (znorm_windows != nullptr && len1 == w_size && len2 == w_size &&
      start1 >= 0 && start2 >= 0 &&
      start1 < (int)znorm_windows->size() &&
      start2 < (int)znorm_windows->size()) {
    return znormed_rms_dist_precomputed(
        znorm_windows->at(start1), znorm_windows->at(start2));
  }

  static thread_local std::vector<double> paa_buf;

  if (len1 == len2) {
    return znormed_rms_dist_equal(series, start1, end1, start2, end2, n_threshold);
  }
  if (len1 < len2) {
    _paa2_range(series, start2, end2, len1, paa_buf);
    return znormed_rms_dist_slice_vec(series, start1, end1, paa_buf, n_threshold);
  }
  _paa2_range(series, start1, end1, len2, paa_buf);
  return znormed_rms_dist_vec_slice(paa_buf, series, start2, end2, n_threshold);
}

double _shrinked_distance(int start1, int end1, int start2, int end2, std::vector<double> *series){

  double res = 0;
  int count = 0;
  int len1 = end1 - start1;
  int len2 = end2 - start2;

  if(len1 == len2) {
    for(int i=0; i<len1; i++){
      res = res + (series->at(start1+i) - series->at(start2+i)) *
        (series->at(start1+i) - series->at(start2+i));
      count++;
    }
    return sqrt(res) / (double) count;
  } else {
    int min_length = std::min(len1, len2);
    for(int i=0; i<min_length; i++){
      res = res + (series->at(start1+i) - series->at(start2+i)) *
        (series->at(start1+i) - series->at(start2+i));
      count++;
    }
    return sqrt(res) / (double) count;
  }
}

double _mean(std::vector<int> *ts, int *start, int *end){
  int sum = 0;
  for(int i=*start; i<*end; i++){
    sum = sum + ts->at(i);
  }
  return (double) sum / (double) (*end - *start);
}

rra_discord_record find_best_rra_discord(std::vector<double> *ts, int w_size,
      std::unordered_map<int, rule_record> *grammar, std::vector<int> *indexes,
      std::vector<rule_interval> *intervals,
      std::unordered_set<int> *global_visited_positions, double n_threshold,
      const std::vector<std::vector<double>>* znorm_windows,
      int seed = -1){

  // *****
  // std::chrono::time_point<std::chrono::system_clock> tstart0, tstart, tend;
  // tstart0 = std::chrono::system_clock::now();
  int distance_calls_counter = 0;

  std::vector<int> visit_array(ts->size(), -1);

  //   for(auto it=intervals.begin(); it!=intervals.end(); ++it){
  //     Rcout << "R" << it->rule_id << " covr " << it->cover << std::endl;
  //   }

  // init variables
  int bestSoFarPosition = -1;
  int bestSoFarLength = -1;
  int bestSoFarRule = 0;
  double bestSoFarDistance = -1;

  // outer loop over all intervals
  for(int i = 0; i < (int) intervals->size(); i++){

    // ****
    // tstart = std::chrono::system_clock::now();

    rule_interval c_interval = intervals->at(i);

    // Rcout << c_interval.rule_id << ", " << c_interval.cover << std::endl;

    auto find = global_visited_positions->find(c_interval.start);
    if(find != global_visited_positions->end()){
      continue;
    }

    // mark the interval location
    std::unordered_set<int> visited_locations;
    visited_locations.reserve(ts->size());

    int markStart = c_interval.start - (c_interval.end - c_interval.start);
    if (markStart < 0) {
      markStart = 0;
    }
    int markEnd = c_interval.end;
    if (markEnd > (int) ts->size()) {
      markEnd = ts->size();
    }
    for(int j=markStart;j<markEnd;j++){
      visited_locations.emplace(j);
    }

    // initialize the distance
    double nn_distance = std::numeric_limits<double>::max();

    // by default, we engage the random search
    bool do_random_search = true;

     //Rcout << " considering interval " << c_interval.start << "-" << c_interval.end <<
    //   " for rule " << c_interval.rule_id <<
    //    ", best so far dist " << bestSoFarDistance << std::endl;

    auto this_rule_occurrences = grammar->at(c_interval.rule_id).rule_intervals;
//     Rcout << "   going to iterate over " << this_rule_occurrences.size() <<
//      " rule occurrences first " << std::endl;

    for(auto it=this_rule_occurrences.begin(); it !=this_rule_occurrences.end(); ++it) {
      // Rcout << "0.0\n";
      int start = indexes->at(it->first);
      // Rcout << "0.1" << start << "\n";
      auto found = visited_locations.find(start);
      if (found == visited_locations.end()) {
        visited_locations.emplace(start);
        int end = indexes->at(it->second) + w_size;
        // Rcout << "    examining a candidate at " << start << "-" <<
         // end << std::endl;
        distance_calls_counter++;
        double dist = _normalized_distance(c_interval.start, c_interval.end,
                                          start, end, *ts, n_threshold,
                                          w_size, znorm_windows);
        // keep track of best so far distance
        if (dist < nn_distance) {
          // Rcout << "    better nn distance found " << dist << std::endl;
          nn_distance = dist;
        }
        if (dist < bestSoFarDistance) {
//           Rcout << "   R " << c_interval.rule_id << ", dist " << dist <<
//            " is less than best so far, breaking off the search" << std::endl;
          do_random_search = false;
          break;
        }
      } else {
        continue;
      }
    }

    // ****
    // tend = std::chrono::system_clock::now();
    // std::chrono::duration<double> elapsed_seconds = tend - tstart;
    // std::cout << "    . pre-search part done in " << elapsed_seconds.count() << "s\n";
    // tstart = std::chrono::system_clock::now();

    if(do_random_search){
//       Rcout << " starting the random search ..." <<
//         " nn dist " << nn_distance << std::endl;
      // Rcout << "visited locations ";
      // for(auto it=visited_locations.begin(); it != visited_locations.end(); ++it){
      //  Rcout << *it << ", ";
      //}
      //Rcout << std::endl;

      // init the visit array
      int cIndex = 0;
      for (unsigned j = 0; j < intervals->size(); j++) {
        rule_interval interval = intervals->at(j);
        auto found = visited_locations.find(interval.start);
        if (found == visited_locations.end()) {
          visit_array[cIndex] = j;
          cIndex++;
          //} else {
          // Rcout << "    - skipped " << interval.start << std::endl;
          //}
        }
      }
      cIndex--;

      // shuffle the visit array
      // only shuffle the filled portion of visit_array
      if (cIndex >= 0) {
        // seed < 0 keeps the historical default-constructed engine (a fixed
        // default seed -- deterministic but not caller-controllable); seed >= 0
        // makes the phase-2 visit order caller-reproducible. Either way the
        // discord result is order-independent; only the distance-call count and
        // search trajectory depend on the order.
        std::default_random_engine rng = (seed < 0)
            ? std::default_random_engine{}
            : std::default_random_engine{(unsigned) seed};
        // std::shuffle(std::begin(visit_array), std::end(visit_array), rng); # old bad code, throws us into -1 zone
        std::shuffle(visit_array.begin(), visit_array.begin() + cIndex + 1, rng);
      }
//      for (int j = cIndex; j > 0; j--) {
//        int index = armaRand() % (j + 1);
//        int a = visit_array[index];
//        visit_array[index] = visit_array[j];
//        visit_array[j] = a;
//      }

      // while there are unvisited locations
      while (cIndex >= 0) {
        int idx = visit_array[cIndex];
        cIndex--;

        // safety check (optional)
        if (idx < 0 || idx >= intervals->size()) continue;

        rule_interval randomInterval = intervals->at(idx);

        distance_calls_counter++;
        double dist = _normalized_distance(
          c_interval.start, c_interval.end,
          randomInterval.start, randomInterval.end, *ts, n_threshold,
          w_size, znorm_windows
        );

        if (dist < bestSoFarDistance) {
          nn_distance = dist;
          break;
        }
        if (dist < nn_distance) {
          nn_distance = dist;
        }
      }

      // while there are unvisited locations
      //while (cIndex >= 0) {

        //rule_interval randomInterval = intervals->at(visit_array[cIndex]);
        //cIndex--;

        //Rcout << "    random candidate " << randomInterval.start << "-" <<
        //  randomInterval.end << ", cindex " << cIndex << std::endl;

        //distance_calls_counter++;
        //double dist = _normalized_distance(c_interval.start, c_interval.end,
        //                randomInterval.start, randomInterval.end, ts);
        //if (dist < bestSoFarDistance) {
          // Rcout << "    dist " << dist <<
          //   " is less than best so far, breaking off the search" << std::endl;
        //  nn_distance = dist;
        //  break;
        //}
        //if (dist < nn_distance) {
          // Rcout << "    better nn distance found " << dist << std::endl;
        //  nn_distance = dist;
        //}

      //}

      // ****
      // tend = std::chrono::system_clock::now();
      // std::chrono::duration<double> elapsed_seconds = tend - tstart;
      // std::cout << "    . random search done in " << elapsed_seconds.count() << "s\n";
      // tstart = std::chrono::system_clock::now();

    } // random search

    if(nn_distance > bestSoFarDistance){
      bestSoFarDistance = nn_distance;
      bestSoFarPosition = c_interval.start;
      bestSoFarLength = c_interval.end - c_interval.start;
      bestSoFarRule = c_interval.rule_id;
      // Rcout << "    updating the discord " << nn_distance << " at " << bestSoFarPosition <<
       // " of length " << bestSoFarLength << " for rule " << bestSoFarRule << std::endl;
    }


  }

  // Rcout << "  RRA, distance calls: " << distance_calls_counter << std::endl;

  rra_discord_record res;
  res.rule = bestSoFarRule;
  res.start = bestSoFarPosition;
  res.end = bestSoFarPosition + bestSoFarLength;
  res.nn_distance = bestSoFarDistance;
  res.distance_calls = distance_calls_counter;

  return res;
}

//' Finds a discord with RRA (Rare Rule Anomaly) algorithm.
//' Usually works the best with higher than that for HOT-SAX sizes of discretization parameters
//' (i.e., PAA and Alphabet sizes).
//'
//' @param series the input timeseries.
//' @param w_size the sliding window size.
//' @param paa_size the PAA size.
//' @param a_size the alphabet size.
//' @param nr_strategy the numerosity reduction strategy ("none", "exact", "mindist").
//' @param n_threshold the normalization threshold.
//' @param discords_num the number of discords to report.
//' @param seed the random seed for the phase-2 search-order shuffle. The default
//' (a negative value) keeps the historical fixed-default engine, which is itself
//' deterministic but not caller-controllable; a non-negative value lets the
//' caller choose the shuffle, so a different reproducible search trajectory (and
//' distance-call count) can be obtained. The reported discords are identical
//' either way -- only the search trajectory depends on it.
//' @useDynLib jmotif
//' @export
//' @references Senin Pavel and Malinchik Sergey,
//' SAX-VSM: Interpretable Time Series Classification Using SAX and Vector Space Model.,
//' Data Mining (ICDM), 2013 IEEE 13th International Conference on.
//' @examples
//' discords = find_discords_rra(ecg0606, 100, 4, 4, "none", 0.01, 1)
//' plot(ecg0606, type = "l", col = "cornflowerblue", main = "ECG 0606")
//' lines(x=c(discords[1,2]:(discords[1,2]+100)),
//'    y=ecg0606[discords[1,2]:(discords[1,2]+100)], col="red")
// [[Rcpp::export]]
Rcpp::DataFrame find_discords_rra(NumericVector series, int w_size, int paa_size,
  int a_size, CharacterVector nr_strategy, double n_threshold = 0.01,
  int discords_num = 3, int seed = -1){

  // *****
  // std::chrono::time_point<std::chrono::system_clock> tstart0, tstart, tend;
  // tstart0 = std::chrono::system_clock::now();
  // tstart = std::chrono::system_clock::now();

  std::vector<double> ts = Rcpp::as<std::vector<double>>(series);
  std::vector<int> visit_array(ts.size(), -1);

  std::unordered_map<int, std::string> sax_map = _sax_via_window(
    ts, w_size, paa_size, a_size, Rcpp::as<std::string>(nr_strategy), n_threshold);

  // ****
  // tend = std::chrono::system_clock::now();
  // std::chrono::duration<double> elapsed_seconds = tend - tstart;
  // Rcout << "  sax conversion: " << elapsed_seconds.count() << "s\n";
  // tstart = std::chrono::system_clock::now();

  // sax_map maps time-series positions to corresponding SAX words
  // to compose the string we need to order keys
  //
  std::vector<int> indexes(sax_map.size());
  int i=0;
  for(auto it = sax_map.begin(); it != sax_map.end(); ++it) {
    indexes[i] = it->first;
    i++;
  }
  sort( indexes.begin(), indexes.end() );

  // Rcout << "  there are " << indexes.size() << " SAX words..." << std::endl;

  // now compose the string
  //
  std::string sax_str;
  for(unsigned i=0; i<indexes.size(); i++){
    sax_str.append(sax_map[ indexes[i] ]).append(" ");
  }
  sax_str.erase( sax_str.end()-1 );

  // ****
  // tend = std::chrono::system_clock::now();
  // elapsed_seconds = tend - tstart;
  // Rcout << "  sax string composition: " << elapsed_seconds.count() << "s\n";

  // grammar
  //
  // *****
  // tstart = std::chrono::system_clock::now();
  std::unordered_map<int, rule_record> grammar = _str_to_repair_grammar(sax_str);

  // ****
  // tend = std::chrono::system_clock::now();
  // elapsed_seconds = tend - tstart;
  // Rcout << "  grammar inferred in: " << elapsed_seconds.count() << "s\n";
  // Rcout << "  there are " << grammar.size() << " RePair rules including R0..." << std::endl;

  // *****
  // tstart = std::chrono::system_clock::now();

  // making intervals and ranking by the rule use
  // meanwhile build the coverage curve
  //
  std::vector<rule_interval> intervals;
  std::vector<int> coverage_array(ts.size(), 0);
  for(auto it = grammar.begin(); it != grammar.end(); ++it) {
    if(0 == it->first){
      continue; // skip R0
    }
    for(auto rit = it->second.rule_intervals.begin(); rit != it->second.rule_intervals.end(); ++rit) {
      int t_start = rit->first;
      int t_end = rit->second;
      // start and end here is for the string tokens, not for time series points
      //
      int start = indexes[t_start];
      int end = indexes[t_end] + w_size;
      // Rcout << start << "_" << end << std::endl;
      // Rcout << " * rule interval " << start << " " << end << std::endl;
      for(int i=start; i<end; ++i){ // rule coverage
        coverage_array[i] = coverage_array[i] + 1;
      }
      rule_interval rr;
      rr.rule_id = it->first;
      rr.start = start;
      rr.end = end;
      // Sort key set at construction (mirrors GrammarViz Java
      // RuleInterval.setCoverage(rule.getRuleIntervals().size()),
      // GrammarVizAnomaly.java:740): RULE FREQUENCY -- the number of
      // occurrences of this interval's grammar rule, NOT per-point coverage.
      // A rule seen once is, by construction, the most anomalous pattern.
      rr.cover = (double) it->second.rule_intervals.size();
      intervals.push_back(rr);
    }
  }

  // we need to examine the coverage curve for continous zero intervals and mark those
  //
  rule_record rec_zero_cover;
  rec_zero_cover.rule_id = -1;
  rec_zero_cover.rule_string = "xxx";
  rec_zero_cover.expanded_rule_string = "xxx";
  bool need_placement = false;
  int start = -1;
  bool in_interval = false;

  // Fix-1 (cross-impl port from the saxpy sibling, commit 822503b): SKIP
  // DEGENERATE zero-coverage runs shorter than a single PAA segment.
  //
  // WHAT: only emit a synthetic rule_id == -1 candidate when the uncovered
  // stretch is at least min_uncovered = max(2, paa_size) points long. A run of
  // length 1 (or any sub-paa_size span) is shorter than one PAA segment, so it
  // is not a discoverable variable-length subsequence -- there is nothing for
  // the RRA NN search to legitimately rank as a discord there. A common trigger
  // is position 0 when RePair leaves the first SAX word as a bare terminal in
  // R0, leaving a single uncovered point.
  //
  // WHY: (a) sub-PAA-segment non-discoverability (above), and (b) cross-impl
  // CONSISTENCY with the saxpy port, which already skips these (saxpy/rra.py,
  // min_uncovered = max(2, paa_size)). In saxpy the run is also numerically
  // pathological: _fast_znorm mean-centers a length-1 span to [0.0], making its
  // distance 0 to *everything*, which both fabricates a spurious nn=0 "discord"
  // and poisons the early-abandon best-so-far, inflating distance calls ~20x
  // (ECG308 100,4,4: saxpy 479,346 calls with the fix). This C++ port does NOT
  // share that catastrophe: _znorm (znorm.cpp) returns the RAW UNCENTERED clone
  // when stdev < threshold, so a length-1 span yields a NONZERO distance, not 0.
  // On the ECG308 100,4,4 row we benchmarked the degenerate span was NOT the
  // winning discord -- but it CAN win on other rows: an A/B build comparison
  // found the C++ original (no Fix-1) reporting rule_id == -1, [0,1) (length 1)
  // as its #1 discord at ecg 100,1,2 and at a seed-99 random walk 100,4,4
  // (nonzero nn ~0.22, non-catastrophic but still a meaningless 1-point span).
  // So Fix-1 does change the result on those rows -- in the INTENDED direction,
  // removing the meaningless top "discord" and surfacing the real ones. Where it
  // was not already winning, this just removes a meaningless extra candidate; in
  // neither case is it fixing a poisoned search. The perf delta is therefore
  // expected to be modest (possibly neutral on rows with no short run) and is to
  // be MEASURED, not assumed.
  //
  // This DELIBERATELY DIVERGES from the authoritative Java GrammarViz
  // (GrammarVizAnomaly.getZeroIntervals, lines 922-939), which keeps EVERY zero
  // run including length 1 (`new RuleInterval(id, start, i, 0)` with no length
  // guard). The divergence is justified by the two points above; the per-point
  // behavior is otherwise identical to Java (cover = 0 still sorts these first).
  //
  // Signedness: the loop counter i is unsigned and start is int, so compute the
  // run length as a signed int (run_len = (int)i - start) before comparing, to
  // avoid an unsigned-promotion comparison.
  int min_uncovered = std::max(2, paa_size);

  for (unsigned i = 0; i < coverage_array.size(); i++) {
    if (0 == coverage_array[i] && !in_interval) {
      start = i;
      in_interval = true;
    }
    if (coverage_array[i] > 0 && in_interval) {

      int run_len = (int) i - start;
      if (run_len >= min_uncovered) {

        need_placement = true;

        rule_interval ri;
        // Zero-coverage stretch: no grammar rule covers it. cover = 0 makes these
        // sort FIRST (rarest-first) -- an uncovered span is maximally anomalous.
        // This matches the authoritative Java GrammarViz exactly: getZeroIntervals
        // builds them as `new RuleInterval(id, start, i, 0)` (coverage == 0,
        // GrammarVizAnomaly.java:933). NOTE this INTENTIONALLY diverges from the
        // saxpy port and from this file's previous post-pass, which both ranked
        // zero spans by their OWN occurrence count (len(zero_rule.intervals));
        // cover = 0 pushes them ahead of every count>=1 rule instead.
        ri.cover = 0;
        ri.start = start;
        ri.end=i;
        ri.rule_id=-1;

        rec_zero_cover.rule_occurrences.push_back(start);
        rec_zero_cover.rule_intervals.push_back(std::make_pair(start, i));

        // NOTE: this interval was previously push_back'd TWICE (once before and
        // once after the rec_zero_cover bookkeeping), duplicating every
        // zero-coverage candidate in the search list -- inflating the candidate
        // set and wasting distance calls on an identical interval. Insert once.
        intervals.push_back(ri);

        // Rcout << " zero coverage from " << start << " to " << i << std::endl;
      }

      // The run always closes when coverage turns positive, whether or not it
      // qualified above -- mirrors saxpy, where `in_interval = False` sits
      // OUTSIDE the `if i - start >= min_uncovered:` length guard.
      in_interval = false;
    }
  }

  if(need_placement) {
    grammar.emplace(-1, rec_zero_cover);
  }

  // Rank the candidate intervals "rarest first" so the early-abandoning NN
  // search builds a strong best-so-far distance quickly.
  //
  // Sort key: RULE FREQUENCY -- the number of occurrences of the interval's
  // grammar rule (a rule seen once is, by construction, the most anomalous
  // pattern). This matches the canonical jMotif GrammarViz (Java) RRA and the
  // saxpy port. We reuse the `cover` field (and the existing sort_intervals
  // comparator, left.cover < right.cover) to carry the frequency, so a small
  // count sorts first. The frequency is now set on each interval at
  // construction time (see the build loop above and the zero-coverage block),
  // so no separate pass is needed here -- the comparator just reads the
  // pre-set field, exactly as the Java reference does
  // (RRAImplementation.findBestDiscordForIntervals sorts on getCoverage()).
  //
  // NOTE -- alternative ordering: MEAN PER-POINT COVERAGE (how many rules
  // overlap each point of the span), jmotif-R's historical key. It is a local
  // density measure rather than a pattern-rarity one, and the two are nearly
  // uncorrelated. To use coverage instead, set cover from the coverage curve
  // with a pass like:
  //     for(auto it=intervals.begin(); it!=intervals.end(); ++it)
  //         it->cover = _mean(&coverage_array, &it->start, &it->end);
  //
  // Tie handling: the key is a small INTEGER frequency, so large tie-groups
  // form (all freq-1 rules, all freq-2 rules, ...). We use std::STABLE_sort to
  // match the two reference impls -- Java's Collections.sort (stable TimSort,
  // single Double.compare key, NO secondary key) and saxpy's list.sort (also
  // stable). stable_sort preserves construction order within a tie-group and
  // removes the run-to-run / STL-version reordering that the previous std::sort
  // could introduce. Note FULL byte-for-byte tie order still cannot match Java:
  // the construction order itself differs because `grammar` is a
  // std::unordered_map (implementation-defined iteration order) whereas Java
  // iterates rules in deterministic id order. We deliberately keep the single
  // frequency key with no secondary tie-break, mirroring all three sibling
  // implementations rather than inventing one here.

  // sort the intervals rare < frequent (stable: preserve construction order on ties)
  std::stable_sort(intervals.begin(), intervals.end(), sort_intervals());


  // ******
  // tend = std::chrono::system_clock::now();
  // elapsed_seconds = tend - tstart;
  // Rcout << "  there are " << intervals.size() <<" rule intervals inferred in "
  //       << elapsed_seconds.count() << std::endl;
  // Rcout << "  top coverage for interval of rule " << intervals[0].rule_id << " starting at "
  //       << intervals[0].start << " ending at " << intervals[0].end << " : " << intervals[0].cover
  //       << std::endl;
  // int last_idx = intervals.size() - 1;
  // Rcout << "  bottom coverage for interval of rule " << intervals[last_idx].rule_id
  //       << " starting at " << intervals[last_idx].start << " ending at " << intervals[last_idx].end
  //       << " : " << intervals[last_idx].cover << std::endl;
  // tstart = std::chrono::system_clock::now();

  // from here on we'll be calling find best discord...
  std::unordered_set<int> global_visited_positions;
  // global_visited_positions.reserve(ts.size());

  int n_windows = (int)ts.size() - w_size + 1;
  std::vector<std::vector<double>> znorm_windows(n_windows);
  for (int wi = 0; wi < n_windows; wi++) {
    _znorm_slice(ts, wi, wi + w_size, n_threshold, znorm_windows[wi]);
  }

  std::vector<rra_discord_record> discords;

  while((int) discords.size() < discords_num){

    // tstart = std::chrono::system_clock::now();

    rra_discord_record d = find_best_rra_discord(&ts, w_size, &grammar,
                              &indexes, &intervals, &global_visited_positions, n_threshold,
                              &znorm_windows, seed);
    // Rcout << d.nn_distance;

    if(d.nn_distance<0){
      break;
    }
    discords.push_back(d);

    // mark the location
    int markStart = d.start - (d.end - d.start);
    if (markStart < 0) {
      markStart = 0;
    }
    int markEnd = d.end;
    if (markEnd > (int) ts.size()) {
      markEnd = ts.size();
    }
    for(int j=markStart;j<markEnd;j++){
      global_visited_positions.emplace(j);
    }

    // *****
    // tend = std::chrono::system_clock::now();
    // elapsed_seconds = tend - tstart;
    // Rcout << "  search for discord " << discords.size() - 1 <<" finished in "
    //       << elapsed_seconds.count() << std::endl;
  }

//   tend = std::chrono::system_clock::now();
//   elapsed_seconds = tend - tstart;
//   std::cout << "discord found in: " << elapsed_seconds.count() << "s\n";

  // make results
  std::vector<int> rule_ids;
  std::vector<int> starts;
  std::vector<int> ends;
  std::vector<int> lengths;
  std::vector<double > nn_distances;
  std::vector<int> distance_calls;


  for(auto it = discords.begin(); it != discords.end(); it++) {
    rule_ids.push_back(it->rule);
    starts.push_back(it->start);
    ends.push_back(it->end);
    lengths.push_back(it->end - it->start);
    nn_distances.push_back(it->nn_distance);
    distance_calls.push_back(it->distance_calls);
  }

  // make results
  return Rcpp::DataFrame::create(
    Named("rule_id") = rule_ids,
    Named("start") = starts,
    Named("end") = ends,
    Named("length") = lengths,
    Named("nn_distance") = nn_distances,
    Named("distance_calls") = distance_calls
  );

}
