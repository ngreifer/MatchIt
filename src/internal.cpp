#include "internal.h"
using namespace Rcpp;

// Rcpp internal functions

//C implementation of tabulate(). Faster than base::tabulate(), but real
//use is in subclass2mmC().

IntegerVector tabulateC_(const IntegerVector& bins,
                         int nbins) {
  int max_bin;

  if (nbins > 0) max_bin = nbins;
  else max_bin = max(na_omit(bins));

  IntegerVector counts(max_bin);
  R_xlen_t n = bins.size();
  for (R_xlen_t i = 0; i < n; i++) {
    if (bins[i] > 0 && bins[i] <= max_bin) {
      counts[bins[i] - 1]++;
    }
  }

  return counts;
}

//Rcpp port of base::which

IntegerVector which(const LogicalVector& x) {
  IntegerVector ind = Range(0, x.size() - 1);
  return ind[x];
}

//Position of `focal` among the sorted unique treatment values, which is how the
//treatment is recoded to 0..g-1 everywhere it is used as an index. `focal` is always
//one of those values, so the `g` returned when it is absent is unreachable.

int recode_focal(int focal,
                 const IntegerVector& unique_treat) {
  int g = unique_treat.size();

  for (int gi = 0; gi < g; gi++) {
    if (unique_treat[gi] == focal) {
      return gi;
    }
  }

  return g;
}

bool antiexact_okay(int aenc,
                    int i,
                    int j,
                    const IntegerMatrix& antiexact_covs) {
  if (aenc == 0) {
    return true;
  }

  //Indexed directly rather than through `.row()`, which would allocate a vector
  //for each of the two rows on every call
  for (int k = 0; k < aenc; k++) {
    if (antiexact_covs(i, k) == antiexact_covs(j, k)) {
      return false;
    }
  }

  return true;
}

bool caliper_covs_okay(int ncc,
                       int i,
                       int j,
                       const NumericMatrix& caliper_covs_mat,
                       const NumericVector& caliper_covs) {
  if (ncc == 0) {
    return true;
  }

  for (int k = 0; k < ncc; k++) {
    if (caliper_covs[k] >= 0) {
      if (std::abs(caliper_covs_mat(i, k) - caliper_covs_mat(j, k)) > caliper_covs[k]) {
        return false;
      }
    }
    else {
      if (std::abs(caliper_covs_mat(i, k) - caliper_covs_mat(j, k)) <= -caliper_covs[k]) {
        return false;
      }
    }
  }

  return true;
}

//Only used in this file
static bool caliper_covs_okay2(int ncc,
                               const NumericVector& cc_ti,
                               int j,
                               const NumericMatrix& caliper_covs_mat,
                               const NumericVector& caliper_covs) {
  if (ncc == 0) {
    return true;
  }

  for (int k = 0; k < ncc; k++) {
    if (caliper_covs[k] >= 0) {
      if (std::abs(cc_ti[k] - caliper_covs_mat(j, k)) > caliper_covs[k]) {
        return false;
      }
    }
    else {
      if (std::abs(cc_ti[k] - caliper_covs_mat(j, k)) <= -caliper_covs[k]) {
        return false;
      }
    }
  }

  return true;
}

bool caliper_dist_okay(bool use_caliper_dist,
                       int i,
                       int j,
                       const NumericVector& distance,
                       double caliper_dist) {
  if (!use_caliper_dist) {
    return true;
  }

  if (caliper_dist >= 0) {
    return std::abs(distance[i] - distance[j]) <= caliper_dist;
  }
  else {
    return std::abs(distance[i] - distance[j]) > -caliper_dist;
  }
}

bool mm_okay(int r,
             int i,
             const IntegerVector& mm_rowi) {

  if (r > 1) {
    //NA_INTEGER is INT_MIN, which never equals a unit index, so the NAs need no
    //separate handling and `na_omit()`, which allocates, is not needed
    R_xlen_t nmm = mm_rowi.size();

    for (R_xlen_t k = 0; k < nmm; k++) {
      if (mm_rowi[k] == i) {
        return false;
      }
    }
  }

  return true;
}

bool exact_okay(bool use_exact,
                int i,
                int j,
                const IntegerVector& exact) {

  if (!use_exact) {
    return true;
  }

  return exact[i] == exact[j];
}

//Squared Euclidean distance between two rows of a matrix. Takes row indices rather
//than vectors because `.row()` allocates, and this runs once per candidate control.
double euc_dist_sq(const NumericMatrix& x,
                   int i,
                   int j) {
  double out = 0;
  int p = x.ncol();

  for (int k = 0; k < p; k++) {
    double tmp = x(i, k) - x(j, k);
    out += tmp * tmp;
  }

  return out;
}

//Return the `ratio` ids with the smallest distances, in increasing order of
//distance. Shared by find_control_vec() and find_control_mat(), which differ in how
//they collect candidates but not in how they choose among them.
std::vector<int> take_closest(std::vector<int> ids,
                              const std::vector<double>& dists,
                              int ratio) {

  int n = ids.size();

  if (n <= 1) {
    return ids;
  }

  if (n <= ratio && std::is_sorted(dists.begin(), dists.end())) {
    return ids;
  }

  std::vector<int> ind(n);
  std::iota(ind.begin(), ind.end(), 0);

  auto dist_less = [&dists](int a, int b) {
    return dists[a] < dists[b];
  };

  //`partial_sort()` and `sort()` order ties differently, so which one runs has to
  //stay as it was for the output to be unchanged
  if (n > ratio) {
    std::partial_sort(ind.begin(), ind.begin() + ratio, ind.end(), dist_less);
    ind.resize(ratio);
  }
  else {
    std::sort(ind.begin(), ind.end(), dist_less);
  }

  std::vector<int> matches_out;
  matches_out.reserve(ind.size());

  for (int i : ind) {
    matches_out.push_back(ids[i]);
  }

  return matches_out;
}

//A counting sort on the stratum codes. Filling each stratum in the order of `units`
//keeps its units in the same relative order they have there, so this is a stable
//sort of `units` by stratum. Units with a missing code, which `method = "cem"` passes
//for units outside every stratum, form one stratum of their own, just as
//`exact_okay()` treats two missing codes as equal.
ExactOrder make_exact_order(const IntegerVector& exact,
                            const IntegerVector& units) {
  ExactOrder out;

  R_xlen_t n = units.size();
  R_xlen_t i;

  int code_min = 0;
  int code_max = -1;
  bool any_code = false;

  for (i = 0; i < n; i++) {
    int code = exact[units[i]];

    if (code == NA_INTEGER) {
      continue;
    }

    if (!any_code) {
      code_min = code_max = code;
      any_code = true;
    }
    else if (code < code_min) {
      code_min = code;
    }
    else if (code > code_max) {
      code_max = code;
    }
  }

  out.code_min = code_min;
  out.na_stratum = code_max - code_min + 1;

  int n_strata = out.na_stratum + 1;

  //`start[e]` becomes the first position of stratum `e`, and `start[n_strata]` is `n`
  std::vector<int> start(n_strata + 1, 0);
  for (i = 0; i < n; i++) {
    start[out.stratum(exact[units[i]]) + 1]++;
  }

  for (int e = 0; e < n_strata; e++) {
    start[e + 1] += start[e];
  }

  out.ord = IntegerVector(n);
  out.pos = IntegerVector(exact.size(), NA_INTEGER);
  out.first = IntegerVector(n_strata);
  out.last = IntegerVector(n_strata);

  for (int e = 0; e < n_strata; e++) {
    out.first[e] = start[e];
    out.last[e] = start[e + 1] - 1;
  }

  std::vector<int> next(start.begin(), start.end() - 1);

  for (i = 0; i < n; i++) {
    int u = units[i];
    int q = next[out.stratum(exact[u])]++;

    out.ord[q] = u;
    out.pos[u] = q;
  }

  return out;
}

std::pair<int, double> find_match_var(const NumericMatrix& mah_covs,
                                      const NumericMatrix& caliper_covs_mat,
                                      const NumericVector& caliper_covs,
                                      const IntegerVector& rows) {
  const int n_mah_covs = mah_covs.ncol();
  const int ncc = caliper_covs_mat.ncol();
  const bool all_rows = rows.size() == 0;

  for (int mci = 0; mci < n_mah_covs; mci++) {
    for (int cci = 0; cci < ncc; cci++) {
      if (caliper_covs[cci] < 0) {
        continue;
      }

      double a;

      if (all_rows) {
        a = get_affine_transformation(caliper_covs_mat.column(cci),
                                      mah_covs.column(mci));
      }
      else {
        NumericVector cal_col = caliper_covs_mat.column(cci);
        NumericVector mah_col = mah_covs.column(mci);

        a = get_affine_transformation(cal_col[rows], mah_col[rows]);
      }

      if (std::abs(a) <= 1e-10) {
        continue;
      }

      return std::make_pair(mci, std::abs(a) * caliper_covs[cci]);
    }
  }

  return std::make_pair(0, static_cast<double>(R_PosInf));
}

//Matching variable: when a caliper covariate is an affine transformation of one of
//the `mah_covs` columns, sorting on that column lets the scan in
//`find_control_mahcovs()` stop as soon as the caliper is exceeded. Only the caliper is
//converted to that column's scale, and only for the matching function's own use; the
//caliper is still enforced on its own scale by `caliper_covs_okay()`, and
//`caliper_covs` and `caliper_covs_mat` belong to the caller and are left alone.
std::vector<double> set_up_mahcovs_scan(StrataScan& scan,
                                        NumericVector& match_var,
                                        IntegerVector& ind_d_ord,
                                        IntegerVector& match_d_ord,
                                        const NumericMatrix& mah_covs,
                                        const NumericMatrix& caliper_covs_mat,
                                        const NumericVector& caliper_covs,
                                        const Nullable<IntegerVector>& strata_,
                                        bool local,
                                        Function& o) {
  const R_xlen_t n = mah_covs.nrow();
  R_xlen_t i;

  scan.use = strata_.isNotNull();
  scan.local = scan.use && local;

  if (scan.use) {
    scan.strata = as<IntegerVector>(strata_);
  }

  if (!scan.local) {
    std::pair<int, double> mv = find_match_var(mah_covs, caliper_covs_mat, caliper_covs,
                                               IntegerVector(0));

    match_var = mah_covs.column(mv.first);

    ind_d_ord = o(match_var);
    ind_d_ord = ind_d_ord - 1; //location of each unit after sorting

    //`ind_d_ord` is a permutation, so its order is just its inverse; computing that
    //directly avoids a second call into R
    match_d_ord = IntegerVector(n);
    for (i = 0; i < n; i++) {
      match_d_ord[ind_d_ord[i]] = static_cast<int>(i);
    }

    if (scan.use) {
      scan.order = make_exact_order(scan.strata, ind_d_ord);
    }

    return std::vector<double>(1, mv.second);
  }

  //Each stratum as a separate match of it alone: the matching variable is chosen from
  //the stratum's units, in their original order, and only they are sorted on it
  const IntegerVector all_units = seq(0, n - 1);
  const ExactOrder by_index = make_exact_order(scan.strata, all_units);
  const R_xlen_t n_strata = by_index.first.size();

  std::vector<int> match_var_col(n_strata, 0);
  std::vector<double> match_var_caliper(n_strata, R_PosInf);

  if (caliper_covs_mat.ncol() > 0) {
    for (R_xlen_t e = 0; e < n_strata; e++) {
      if (by_index.first[e] > by_index.last[e]) {
        continue;
      }

      const IntegerVector rows = by_index.ord[Range(by_index.first[e], by_index.last[e])];
      std::pair<int, double> mv = find_match_var(mah_covs, caliper_covs_mat, caliper_covs,
                                                 rows);

      match_var_col[e] = mv.first;
      match_var_caliper[e] = mv.second;
    }
  }

  match_var = NumericVector(n);
  for (i = 0; i < n; i++) {
    match_var[i] = mah_covs(i, match_var_col[by_index.stratum(scan.strata[i])]);
  }

  //Sorting on the stratum and then on its matching variable puts each stratum's units
  //in exactly the order sorting that stratum alone would, because the sort is stable:
  //ties keep the units' original order in both.
  ind_d_ord = o(scan.strata, match_var);
  ind_d_ord = ind_d_ord - 1;

  scan.order = make_exact_order(scan.strata, ind_d_ord);
  match_d_ord = scan.order.pos;

  return match_var_caliper;
}

std::vector<int> find_control_vec(int t_id,
                                  const IntegerVector& ind_d_ord,
                                  const IntegerVector& match_d_ord,
                                  const IntegerVector& treat,
                                  const NumericVector& distance,
                                  const LogicalVector& eligible,
                                  int gi,
                                  int r,
                                  const IntegerVector& mm_rowi_,
                                  int ncc,
                                  const NumericMatrix& caliper_covs_mat,
                                  const NumericVector& caliper_covs,
                                  double caliper_dist,
                                  bool use_exact,
                                  const IntegerVector& exact,
                                  int aenc,
                                  const IntegerMatrix& antiexact_covs,
                                  const IntegerVector& first_control,
                                  const IntegerVector& last_control,
                                  const ExactOrder& exact_order,
                                  int ratio,
                                  int prev_start) {

  //With `exact`, only the treated unit's stratum is scanned, in `exact_order`, which
  //holds each stratum's units in the order they have in `ind_d_ord`. Scanning the
  //whole sample instead steps over every other stratum's units to reach them, which
  //costs time in proportion to the number of strata.
  const IntegerVector& scan_ord = use_exact ? exact_order.ord : ind_d_ord;
  const IntegerVector& scan_pos = use_exact ? exact_order.pos : match_d_ord;

  int ii = scan_pos[t_id];

  IntegerVector mm_rowi;
  std::vector<int> possible_starts;

  if (r > 1) {
    mm_rowi = na_omit(mm_rowi_);
    mm_rowi = mm_rowi[as<IntegerVector>(treat[mm_rowi]) == gi];
    possible_starts.reserve(mm_rowi.size() + 2);

    for (int mmi : mm_rowi) {
      possible_starts.push_back(scan_pos[mmi]);
    }
  }
  else {
    possible_starts.reserve(2);
  }

  if (prev_start >= 0) {
    possible_starts.push_back(scan_pos[prev_start]);
  }

  int iil, iir;
  double min_dist;

  if (possible_starts.empty()) {
    iil = ii;
    iir = ii;
    min_dist = 0;
  }
  else {
    possible_starts.push_back(ii);

    iil = *std::min_element(possible_starts.begin(), possible_starts.end());
    iir = *std::max_element(possible_starts.begin(), possible_starts.end());

    if (iil == ii) {
      min_dist = std::abs(distance[t_id] - distance[scan_ord[iir]]);
    }
    else if (iir == ii) {
      min_dist = std::abs(distance[t_id] - distance[scan_ord[iil]]);
    }
    else {
      min_dist = std::max(std::abs(distance[t_id] - distance[scan_ord[iil]]),
                          std::abs(distance[t_id] - distance[scan_ord[iir]]));
    }
  }

  if (caliper_dist <= 0 && min_dist < -caliper_dist) {
    min_dist = -caliper_dist;
  }

  int min_ii, max_ii;

  if (use_exact) {
    int e = exact_order.stratum(exact[t_id]);
    min_ii = exact_order.first[e];
    max_ii = exact_order.last[e];
  }
  else {
    min_ii = first_control[gi];
    max_ii = last_control[gi];
  }

  //Positions in `ind_d_ord` of the two starting points, for choosing between the
  //sides below when only the stratum is scanned
  const int gl0 = match_d_ord[scan_ord[iil]];
  const int gr0 = match_d_ord[scan_ord[iir]];

  double di = distance[t_id];

  NumericVector cc_ti;
  if (ncc > 0) {
    cc_ti = caliper_covs_mat.row(t_id);
  }

  bool l_stop = false;
  bool r_stop = false;

  double dist_c;

  std::vector<int> potential_matches_id;
  potential_matches_id.reserve(2 * ratio);
  std::vector<double> potential_matches_dist;
  potential_matches_dist.reserve(2 * ratio);

  int num_matches_l = 0;
  int num_matches_r = 0;

  int iz;
  bool left = false;
  int num_closer_than_dist_c;

  while (!l_stop || !r_stop) {
    if (l_stop) {
      left = false;
    }
    else if (r_stop) {
      left = true;
    }
    else if (use_exact) {
      //Taking the sides in the order a scan of the whole sample would reach them keeps
      //the candidates, and the order they are listed in for take_closest(), the same.
      //That scan alternates single steps, so it reaches a unit `a` places left of its
      //left start before one `b` places right of its right start exactly when
      //a <= b. The other strata's units it also steps over can only stop a side, by
      //the caliper or by being farther than `ratio` candidates already found, and both
      //of those also stop the side at the next unit of this stratum, which is farther
      //still and is reached with no fewer candidates found.
      if (iil <= min_ii || num_matches_l == ratio) {
        left = true;
      }
      else if (iir >= max_ii || num_matches_r == ratio) {
        left = false;
      }
      else {
        left = gl0 - match_d_ord[scan_ord[iil - 1]] <= match_d_ord[scan_ord[iir + 1]] - gr0;
      }
    }
    else {
      left = !left;
    }

    if (left) {
      if (iil <= min_ii || num_matches_l == ratio) {
        l_stop = true;
        continue;
      }

      iil -= 1;
      iz = scan_ord[iil];
    }
    else {
      if (iir >= max_ii || num_matches_r == ratio) {
        r_stop = true;
        continue;
      }

      iir += 1;
      iz = scan_ord[iir];
    }

    if (!eligible[iz]) {
      continue;
    }

    if (treat[iz] != gi) {
      continue;
    }

    if (!mm_okay(r, iz, mm_rowi)) {
      continue;
    }

    dist_c = std::abs(di - distance[iz]);

    if (caliper_dist >= 0) {
      if (dist_c > caliper_dist) {
        if (left) {
          l_stop = true;
        }
        else {
          r_stop = true;
        }
        continue;
      }
    }
    else if (dist_c <= -caliper_dist) {
      continue;
    }

    if (dist_c < min_dist) {
      continue;
    }

    //If current dist is worse than ratio dists, continue
    if (potential_matches_id.size() >= static_cast<size_t>(ratio)) {
      num_closer_than_dist_c = 0;
      for (double d : potential_matches_dist) {
        if (d < dist_c) {
          num_closer_than_dist_c++;
          if (num_closer_than_dist_c == ratio) {
            break;
          }
        }
      }

      if (num_closer_than_dist_c >= ratio) {
        if (left) {
          l_stop = true;
        }
        else {
          r_stop = true;
        }
        continue;
      }
    }

    if (!exact_okay(use_exact, t_id, iz, exact)) {
      continue;
    }

    if (!antiexact_okay(aenc, t_id, iz, antiexact_covs)) {
      continue;
    }

    if (!caliper_covs_okay2(ncc, cc_ti, iz, caliper_covs_mat, caliper_covs)) {
      continue;
    }

    potential_matches_id.push_back(iz);
    potential_matches_dist.push_back(dist_c);

    if (left) {
      num_matches_l++;
      if (num_matches_l == ratio) {
        l_stop = true;
      }
    }
    else {
      num_matches_r++;
      if (num_matches_r == ratio) {
        r_stop = true;
      }
    }
  }

  return take_closest(std::move(potential_matches_id), potential_matches_dist, ratio);
}

std::vector<int> find_control_mahcovs(int t_id,
                                      const IntegerVector& ind_d_ord,
                                      const IntegerVector& match_d_ord,
                                      const NumericVector& match_var,
                                      double match_var_caliper,
                                      const IntegerVector& treat,
                                      const NumericVector& distance,
                                      const LogicalVector& eligible,
                                      int gi,
                                      int r,
                                      const IntegerVector& mm_rowi,
                                      const NumericMatrix& mah_covs,
                                      int ncc,
                                      const NumericMatrix& caliper_covs_mat,
                                      const NumericVector& caliper_covs,
                                      bool use_caliper_dist,
                                      double caliper_dist,
                                      bool use_exact,
                                      const IntegerVector& exact,
                                      int aenc,
                                      const IntegerMatrix& antiexact_covs,
                                      const StrataScan& scan,
                                      int ratio) {

  //With strata, only the treated unit's stratum is searched, in `scan.order`
  const IntegerVector& scan_ord = scan.use ? scan.order.ord : ind_d_ord;
  const IntegerVector& scan_pos = scan.use ? scan.order.pos : match_d_ord;

  int ii = scan_pos[t_id];

  int iil, iir;

  iil = ii;
  iir = ii;

  int min_ii, max_ii;

  if (scan.use) {
    int e = scan.order.stratum(scan.strata[t_id]);
    min_ii = scan.order.first[e];
    max_ii = scan.order.last[e];
  }
  else {
    min_ii = 0;
    max_ii = match_d_ord.size() - 1;
  }

  //Position in `ind_d_ord` of the treated unit, where a search of the whole sample
  //starts; used to choose between the sides below
  const bool global_sides = scan.use && !scan.local;
  const int g0 = global_sides ? match_d_ord[t_id] : 0;

  bool l_stop = false;
  bool r_stop = false;

  double dist_c;

  std::vector<std::pair<int, double>> potential_matches;
  potential_matches.reserve(ratio);

  std::pair<int,double> new_match;

  int num_matches_l = 0;
  int num_matches_r = 0;

  double mv_i = match_var[t_id];
  double mv_dist;

  int iz;
  bool left = false;

  auto dist_comp = [](const std::pair<int, double>& a,
                      const std::pair<int, double>& b) {
    return a.second < b.second;
  };

  while (!l_stop || !r_stop) {
    if (l_stop) {
      left = false;
    }
    else if (r_stop) {
      left = true;
    }
    else if (global_sides) {
      //As in find_control_vec(): a search of the whole sample alternates single steps
      //from the treated unit, so it reaches a unit `a` places to its left before one
      //`b` places to its right exactly when a <= b. The other strata's units it also
      //steps over can only stop a side, and the next unit of this stratum would stop
      //it too. A separate match of the stratum (`scan.local`) alternates steps within
      //the stratum instead, which is the plain alternation below.
      if (iil <= min_ii) {
        left = true;
      }
      else if (iir >= max_ii) {
        left = false;
      }
      else {
        left = g0 - match_d_ord[scan_ord[iil - 1]] <= match_d_ord[scan_ord[iir + 1]] - g0;
      }
    }
    else {
      left = !left;
    }

    if (left) {
      if (iil <= min_ii || num_matches_l == ratio) {
        l_stop = true;
        continue;
      }

      iil -= 1;
      iz = scan_ord[iil];
    }
    else {
      if (iir >= max_ii || num_matches_r == ratio) {
        r_stop = true;
        continue;
      }

      iir += 1;
      iz = scan_ord[iir];
    }

    if (!eligible[iz]) {
      continue;
    }

    if (treat[iz] != gi) {
      continue;
    }

    if (!mm_okay(r, iz, mm_rowi)) {
      continue;
    }

    mv_dist = std::abs(mv_i - match_var[iz]);

    if (match_var_caliper >= 0) {
      if (mv_dist > match_var_caliper) {
        if (left) {
          l_stop = true;
        }
        else {
          r_stop = true;
        }
        continue;
      }
    }
    else if (mv_dist <= -match_var_caliper) {
      continue;
    }

    mv_dist = mv_dist * mv_dist;

    //If current dist is worse than ratio dists, continue
    if (potential_matches.size() == static_cast<size_t>(ratio)) {
      if (potential_matches.back().second < mv_dist) {
        if (left) {
          l_stop = true;
        }
        else {
          r_stop = true;
        }

        continue;
      }
    }

    if (!exact_okay(use_exact, t_id, iz, exact)) {
      continue;
    }

    if (!caliper_dist_okay(use_caliper_dist, t_id, iz, distance, caliper_dist)) {
      continue;
    }

    if (!antiexact_okay(aenc, t_id, iz, antiexact_covs)) {
      continue;
    }

    if (!caliper_covs_okay(ncc, t_id, iz, caliper_covs_mat, caliper_covs)) {
      continue;
    }

    dist_c = euc_dist_sq(mah_covs, t_id, iz);

    if (!std::isfinite(dist_c)) {
      continue;
    }

    new_match = std::pair<int,double>(iz, dist_c);

    if (potential_matches.empty()) {
      potential_matches.push_back(new_match);
    }
    else if (dist_c > potential_matches.back().second) {
      if (potential_matches.size() == static_cast<size_t>(ratio)) {
        continue;
      }

      potential_matches.push_back(new_match);
    }
    else if (ratio == 1) {
      potential_matches[0] = new_match;
    }
    else {
      if (potential_matches.size() == static_cast<size_t>(ratio)) {
        potential_matches.pop_back();
      }

      if (dist_c > potential_matches.back().second) {
        potential_matches.push_back(new_match);
      }
      else {
        potential_matches.insert(std::lower_bound(potential_matches.begin(), potential_matches.end(),
                                                  new_match, dist_comp),
                                                  new_match);
      }
    }
  }

  std::vector<int> matches_out;
  matches_out.reserve(potential_matches.size());

  for (const auto& p : potential_matches) {
    matches_out.push_back(p.first);
  }

  return matches_out;
}

std::vector<int> find_control_mat(int t_id,
                                  const IntegerVector& treat,
                                  const IntegerVector& ind_non_focal,
                                  const IntegerVector& ind_match,
                                  const NumericMatrix& distance_mat,
                                  int t_row,
                                  const LogicalVector& eligible,
                                  int gi,
                                  int r,
                                  const IntegerVector& mm_rowi,
                                  int ncc,
                                  const NumericMatrix& caliper_covs_mat,
                                  const NumericVector& caliper_covs,
                                  double caliper_dist,
                                  bool use_exact,
                                  const IntegerVector& exact,
                                  int aenc,
                                  const IntegerMatrix& antiexact_covs,
                                  const StrataScan& scan,
                                  int ratio) {

  int c_id_i;
  double dist_c;

  std::vector<int> potential_matches_id;

  if (ratio < 1) {
    return potential_matches_id;
  }

  std::vector<double> potential_matches_dist;
  double max_dist = R_PosInf;

  //The controls are visited in column order. With strata, only the treated unit's
  //stratum's controls are, which `scan.order` holds in column order; the others would
  //all be passed over, so the candidates are the same either way.
  R_xlen_t k_first = 0;
  R_xlen_t k_last = distance_mat.ncol();

  if (scan.use) {
    int e = scan.order.stratum(scan.strata[t_id]);

    if (e < 0) {
      return potential_matches_id;
    }

    k_first = scan.order.first[e];
    k_last = scan.order.last[e] + 1;
  }

  potential_matches_id.reserve(k_last - k_first);
  potential_matches_dist.reserve(k_last - k_first);

  for (R_xlen_t k = k_first; k < k_last; k++) {
    R_xlen_t c = scan.use ? ind_match[scan.order.ord[k]] : k;

    dist_c = distance_mat(t_row, c);

    if (potential_matches_id.size() >= static_cast<size_t>(ratio)) {
      if (dist_c > max_dist) {
        continue;
      }
    }

    if (caliper_dist >= 0) {
      if (dist_c > caliper_dist) {
        continue;
      }
    }
    else {
      if (dist_c <= -caliper_dist) {
        continue;
      }
    }

    if (!std::isfinite(dist_c)) {
      continue;
    }

    c_id_i = ind_non_focal[c];

    if (!eligible[c_id_i]) {
      continue;
    }

    if (treat[c_id_i] != gi) {
      continue;
    }

    if (!mm_okay(r, c_id_i, mm_rowi)) {
      continue;
    }

    if (!exact_okay(use_exact, t_id, c_id_i, exact)) {
      continue;
    }

    if (!antiexact_okay(aenc, t_id, c_id_i, antiexact_covs)) {
      continue;
    }

    if (!caliper_covs_okay(ncc, t_id, c_id_i, caliper_covs_mat, caliper_covs)) {
      continue;
    }

    potential_matches_id.push_back(c_id_i);
    potential_matches_dist.push_back(dist_c);

    if (potential_matches_id.size() == 1) {
      max_dist = dist_c;
    }
    else if (dist_c > max_dist) {
      max_dist = dist_c;
    }
  }

  return take_closest(std::move(potential_matches_id), potential_matches_dist, ratio);
}

double max_finite(const NumericVector& x) {
  double m = NA_REAL;

  R_xlen_t n = x.size();
  R_xlen_t i;
  bool found = false;

  //Find first finite value
  for (i = 0; i < n; i++) {
    if (std::isfinite(x[i])) {
      m = x[i];
      found = true;
      break;
    }
  }

  //If none found, return NA
  if (!found) {
    return m;
  }

  //Find largest finite value
  for (R_xlen_t j = i + 1; j < n; j++) {
    if (!std::isfinite(x[j])) {
      continue;
    }

    if (x[j] > m) {
      m = x[j];
    }
  }

  return m;
}

double min_finite(const NumericVector& x) {
  double m = NA_REAL;

  R_xlen_t n = x.size();
  R_xlen_t i;
  bool found = false;

  //Find first finite value
  for (i = 0; i < n; i++) {
    if (std::isfinite(x[i])) {
      m = x[i];
      found = true;
      break;
    }
  }

  //If none found, return NA
  if (!found) {
    return m;
  }

  //Find smallest finite value
  for (R_xlen_t j = i + 1; j < n; j++) {
    if (!std::isfinite(x[j])) {
      continue;
    }

    if (x[j] < m) {
      m = x[j];
    }
  }

  return m;
}

void update_first_and_last_control(IntegerVector first_control,
                                   IntegerVector last_control,
                                   const IntegerVector& ind_d_ord,
                                   const LogicalVector& eligible,
                                   const IntegerVector& treat,
                                   int gi) {
  R_xlen_t c;

  // Update first_control
  if (!eligible[ind_d_ord[first_control[gi]]]) {
    for (c = first_control[gi] + 1; c <= last_control[gi]; c++) {
      if (eligible[ind_d_ord[c]]) {
        if (treat[ind_d_ord[c]] == gi) {
          first_control[gi] = c;
          break;
        }
      }
    }
  }

  // Update last_control
  if (!eligible[ind_d_ord[last_control[gi]]]) {
    for (c = last_control[gi] - 1; c >= first_control[gi]; c--) {
      if (eligible[ind_d_ord[c]]) {
        if (treat[ind_d_ord[c]] == gi) {
          last_control[gi] = c;
          break;
        }
      }
    }
  }
}

double get_affine_transformation(const NumericVector& x,
                                 const NumericVector& y,
                                 double tol) {
  R_xlen_t n = x.size();
  R_xlen_t i;

  if (n != y.size() || n < 2) {
    return 0.0; // Need at least two points for a meaningful check
  }

  // Compute means
  double mean_x = mean(x);
  double mean_y = mean(y);

  // Compute a (scaling factor)
  double num = 0.0, denom = 0.0;
  double x_diff, y_diff;
  for (i = 0; i < n; i++) {
    x_diff = x[i] - mean_x;
    y_diff = y[i] - mean_y;

    num += x_diff * y_diff;
    denom += x_diff * x_diff;
  }

  if (std::abs(denom) < tol || std::abs(num) < tol) {
    return 0.0;
  }

  double a = num / denom;
  double b = mean_y - a * mean_x;

  // Verify if y is reconstructed correctly within tolerance
  for (i = 0; i < n; i++) {
    if (std::abs(a * x[i] + b - y[i]) > tol) {
      return 0.0;
    }
  }

  return a;
}
