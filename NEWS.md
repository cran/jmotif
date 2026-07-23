# jmotif 1.3.2

* **Fixed:** RRA phase-2 `visit_array` buffer sized to grammar interval count (was
  series length); could overflow and **segfault** when `seed >= 0` enabled the
  seeded random-search path. Tier-B tests now exercise `seed = 0` on ecg0606.
* **Tests:** tier-B RRA conformance pins on ecg0606 — exact span and NN distance at
  w100 and w120 (`seed = 0`); w150/p7/a4 asserts the top discord is not a
  zero-distance boundary gap.
* **Docs:** README lists CRAN install first; fix header anchors for SAX-VSM and
  RePair; rename RePair walkthrough to §7.1; link saxpy **2.0.1** on PyPI;
  roxygen `\examples{}` on `sax_via_window`, SAX-VSM bag/TF·IDF helpers, and
  `cosine_sim`; clarify `cosine_dist()` title and `find_discords_rra()` defaults
  in `.Rd` files; Codecov badge URL → `app.codecov.io` (CRAN incoming check).
* **Packaging:** ship `NEWS.md` in the tarball; exclude agent audit HTML from
  builds; refresh `DESCRIPTION` summary; roxygen2 8.0.0 metadata.

# jmotif 1.3.1

* Memory and performance (Phases 1–4):
  - Fix memory leaks in RePair grammar (`rule_record` by-value), RRA, and
    `cosine_sim` (Valgrind-clean on long runs).
  - HOT-SAX: precompute z-normed sliding windows (~1.4× on ecg0606, ~17× on
    longer series); identical discord positions and distances.
  - SAX sliding window: `_paa2` in-place fractional PAA, `_znorm_slice`,
    unified `_sax_via_window` path (~2.3× on bag construction).
  - RRA: fused normalized distance, precomputed `w_size` z-norm windows
    (~1.3× on ecg0606); identical results and distance-call counts.
* Safety: `idx_to_letter()` validates index range [1, 26].
* Internal: use `_is_equal_mindist(std::string)` in SAX numerosity reduction;
  remove unused HOT-SAX state; clarify `cosine_dist()` roxygen.

# Version 1.3.0
  Matrix Profile / MASS convention:
  - znorm now uses POPULATION standard deviation (divide by n) instead of
    sample std (n-1). Each normalized window then has empirical variance
    exactly 1, the assumption behind SAX's equiprobable Gaussian breakpoints.
  - SAX symbol assignment maps a value exactly on a breakpoint to the symbol
    ABOVE the cut (cuts[j] <= ts[i]), matching the jmotif Java original and
    saxpy. Previously a strict `<` mapped on-cut values to the symbol below.
  - PAA: fixed a floating-point boundary bug in the fractional PAA (_paa2,
    the public paa() / sax_by_chunking path). The final segment break
    paa_num * (len/paa_num) can round above len in IEEE-754 (e.g. 7*(29/7) =
    29.000000000000004), which dropped the last sample from the final segment
    and corrupted its mean on ~5-6% of non-divisible (n, paa_num) shapes. The
    final break is now snapped to len, matching the saxpy / Java fractional PAA.
  NOTE: all three are behavioral changes -- SAX words for windows whose
  normalized PAA values fall near a breakpoint may differ from 1.2.x output.
* HOT-SAX and brute-force discord discovery now measure the nearest-neighbour
  distance on Z-NORMED windows (the SAX premise of shape similarity):
  - find_discords_hotsax: the three internal _znorm() calls had been left
    commented out, so distances were computed on RAW subsequences -- it reported
    raw-amplitude outliers (e.g. ecg0606 pos 411 / 1.50) instead of true
    z-normed shape discords (pos 430 / 5.28). The windows are now z-normed and
    compared with an early-abandoning Euclidean distance.
  - find_discords_brute_force now z-normalizes too and gained an `n_threshold`
    parameter (default 0.01); previously it had no threshold and ran on raw
    subsequences. It is again a valid z-normed reference for HOT-SAX.
  - Both now use a deterministic lowest-index tie-break and a +/-(w-1) exclusion
    zone, so results are reproducible and match saxpy / Java exactly
    (ecg0606 top discords: 430/5.279, 318/4.176, 2080/2.393).
  NOTE: behavioral change -- discord positions/distances differ from 1.2.x,
  which reported raw-distance discords.

# Version 1.2.1
* Fixed CRAN WARN issue by removing C++ requirement from description

# Version 1.2.0
* Fixed blocking C++ issues with RcppArmaddillo 
* Removed hard requirement for C++ 11
* Fixed R exports into the namespace
* Fixed occasional RRA discord exception when hitting outside of interval
* Fixed broken tests
* Improved some code for efficiency

# Version 1.1.1
* Fixed std::bind2nd deprecation warning.

# Version 1.1.0
* Fixed deprecation and comparison warnings, C++ 11 compatibility.

# Version 1.0.3
* Fixed vector overruns in RRA.

# Version 1.0.2.9000

* Pre-release development version.

* Added the Re-Pair grammatical inference algorithm implementation.

* Added grammatical "rule coverage" curve implementation facilitating 
  a time series anomaly discovery.

* Added RRA (Rare Rule Anomaly) algorithm implementation for grammar-based variable 
  length time series anomaly discovery.

# Version 1.0.2

* Added the HOT-SAX time series discord finding algorithm implementation.

* Added the brute-force time series discord finding algorithm implementation 
  based on the early-abandoning Euclidean distance.

# Version 1.0.1

* CRAN submission bug fixes.

# Version 1.0.0

* Initial submission.
