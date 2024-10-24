// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include "eta_progress_bar.h"
#include <Rcpp.h>
#include "internal.h"
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
IntegerMatrix nn_matchC_vec(const IntegerVector& treat_,
                            const IntegerVector& ord,
                            const IntegerVector& ratio,
                            const LogicalVector& discarded,
                            const int& reuse_max,
                            const int& focal_,
                            const NumericVector& distance,
                            const Nullable<IntegerMatrix>& exact_ = R_NilValue,
                            const Nullable<double>& caliper_dist_ = R_NilValue,
                            const Nullable<NumericVector>& caliper_covs_ = R_NilValue,
                            const Nullable<NumericMatrix>& caliper_covs_mat_ = R_NilValue,
                            const Nullable<IntegerMatrix>& antiexact_covs_ = R_NilValue,
                            const Nullable<IntegerVector>& unit_id_ = R_NilValue,
                            const bool& disl_prog = false) {

  IntegerVector unique_treat = unique(treat_);
  std::sort(unique_treat.begin(), unique_treat.end());
  int g = unique_treat.size();
  IntegerVector treat = match(treat_, unique_treat) - 1;
  int focal;
  for (focal = 0; focal < g; focal++) {
    if (unique_treat[focal] == focal_) {
      break;
    }
  }

  R_xlen_t n = treat.size();
  IntegerVector ind = Range(0, n - 1);

  R_xlen_t i;
  int gi;
  IntegerVector indt(n);
  IntegerVector indt_sep(g + 1);
  IntegerVector indt_tmp;
  IntegerVector nt(g);
  IntegerVector ind_match(n);
  ind_match.fill(NA_INTEGER);

  IntegerVector times_matched(n);
  times_matched.fill(0);
  LogicalVector eligible = !discarded;

  IntegerVector g_c = Range(0, g - 1);
  g_c = g_c[g_c != focal];

  for (gi = 0; gi < g; gi++) {
    nt[gi] = sum(treat == gi);
  }

  int nf = nt[focal];

  indt_sep[0] = 0;

  for (gi = 0; gi < g; gi++) {
    indt_sep[gi + 1] = indt_sep[gi] + nt[gi];

    indt_tmp = ind[treat == gi];

    for (i = 0; i < nt[gi]; i++) {
      indt[indt_sep[gi] + i] = indt_tmp[i];
      ind_match[indt_tmp[i]] = i;
    }
  }

  IntegerVector ind_focal = indt[Range(indt_sep[focal], indt_sep[focal + 1] - 1)];

  IntegerVector times_matched_allowed(n);
  times_matched_allowed.fill(reuse_max);
  times_matched_allowed[ind_focal] = ratio;

  IntegerVector n_eligible(g);
  for (i = 0; i < n; i++) {
    if (eligible[i]) {
      n_eligible[treat[i]]++;
    }
  }

  int max_ratio = max(ratio);

  // Output matrix with sample indices of control units
  IntegerMatrix mm(nf, max_ratio);
  mm.fill(NA_INTEGER);
  CharacterVector lab = treat_.names();

  //Use base::order() because faster than Rcpp implementation of order()
  Function o("order");

  IntegerVector ind_d_ord = o(distance);
  ind_d_ord = ind_d_ord - 1; //location of each unit after sorting

  IntegerVector match_d_ord = match(ind, ind_d_ord) - 1;

  IntegerVector last_control(g);
  last_control.fill(n - 1);
  IntegerVector first_control(g);
  first_control.fill(0);

  //exact
  bool use_exact = false;
  IntegerVector exact;
  if (exact_.isNotNull()) {
    exact = as<IntegerVector>(exact_);
    use_exact = true;
  }

  //caliper_dist
  double caliper_dist;
  if (caliper_dist_.isNotNull()) {
    caliper_dist = as<double>(caliper_dist_);
  }
  else {
    caliper_dist = max_finite(distance) - min_finite(distance) + 1;
  }

  //caliper_covs
  NumericVector caliper_covs;
  NumericMatrix caliper_covs_mat;
  int ncc = 0;
  if (caliper_covs_.isNotNull()) {
    caliper_covs = as<NumericVector>(caliper_covs_);
    caliper_covs_mat = as<NumericMatrix>(caliper_covs_mat_);
    ncc = caliper_covs_mat.ncol();
  }

  //antiexact
  IntegerMatrix antiexact_covs;
  int aenc = 0;
  if (antiexact_covs_.isNotNull()) {
    antiexact_covs = as<IntegerMatrix>(antiexact_covs_);
    aenc = antiexact_covs.ncol();
  }

  //unit_id
  IntegerVector unit_id;
  bool use_unit_id = false;
  if (unit_id_.isNotNull()) {
    unit_id = as<IntegerVector>(unit_id_);
    use_unit_id = true;
  }

  IntegerVector matches_i(g);
  int k_total;

  //progress bar
  int prog_length = sum(ratio) + 1;
  ETAProgressBar pb;
  Progress p(prog_length, disl_prog, pb);

  int r, t_id_t_i, t_id_i, c_id_i, c;
  IntegerVector ck_;
  // bool check = true;

  int counter = -1;

  for (r = 1; r <= max_ratio; r++) {
    for (i = 0; i < nf && max(as<IntegerVector>(n_eligible[g_c])) > 0; i++) {
      // i: generic looping index
      // t_id_t_i; index of treated unit to match among treated units
      // t_id_i: index of treated unit to match among all units
      counter++;
      if (counter % 200 == 0) Rcpp::checkUserInterrupt();

      t_id_t_i = ord[i] - 1;
      t_id_i = ind_focal[t_id_t_i];

      if (r > times_matched_allowed[t_id_i]) {
        continue;
      }

      p.increment();

      if (!eligible[t_id_i]) {
        continue;
      }

      k_total = 0;

      for (int gi : g_c) {
        update_first_and_last_control(first_control,
                                      last_control,
                                      ind_d_ord,
                                      eligible,
                                      treat,
                                      gi);

        c_id_i = find_both(t_id_i,
                           ind_d_ord,
                           match_d_ord,
                           treat,
                           distance,
                           eligible,
                           gi,
                           r,
                           mm.row(t_id_t_i),
                           ncc,
                           caliper_covs_mat,
                           caliper_covs,
                           caliper_dist,
                           use_exact,
                           exact,
                           aenc,
                           antiexact_covs,
                           first_control,
                           last_control);

        if (c_id_i < 0) {
          if (r == 1) {
            k_total = 0;
            break;
          }
          continue;
        }

        matches_i[k_total] = c_id_i;
        k_total++;
      }

      if (k_total == 0) {
        eligible[t_id_i] = false;
        n_eligible[focal]--;
        continue;
      }

      for (c = 0; c < k_total; c++) {
        mm(t_id_t_i, sum(!is_na(mm(t_id_t_i, _)))) = matches_i[c];
      }

      matches_i[k_total] = t_id_i;

      ck_ = matches_i[Range(0, k_total)];

      if (use_unit_id) {
        ck_ = which(!is_na(match(unit_id, as<IntegerVector>(unit_id[ck_]))));
      }

      for (int ck : ck_) {
        if (!eligible[ck]) {
          continue;
        }

        times_matched[ck]++;
        if (times_matched[ck] >= times_matched_allowed[ck]) {
          eligible[ck] = false;
          n_eligible[treat[ck]]--;
        }
      }
    }
  }

  p.update(prog_length);

  mm = mm + 1;
  rownames(mm) = lab[ind_focal];

  return mm;
}
