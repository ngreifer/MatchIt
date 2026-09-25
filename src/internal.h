#ifndef INTERNAL_H
#define INTERNAL_H

#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <utility>
#include <vector>

//No `using namespace Rcpp;` here: this header is included by every translation unit
//in the package, and each of them has its own `using` directive.

Rcpp::IntegerVector tabulateC_(const Rcpp::IntegerVector& bins,
                               int nbins = 0);

Rcpp::IntegerVector which(const Rcpp::LogicalVector& x);

int recode_focal(int focal,
                 const Rcpp::IntegerVector& unique_treat);

//Units grouped by `exact` stratum, each stratum's units in the order they are given
//in, so that a search for controls can visit a treated unit's stratum alone. Built by
//make_exact_order(); empty when there is no `exact`.
struct ExactOrder {
  Rcpp::IntegerVector ord;   //units, stratum by stratum
  Rcpp::IntegerVector pos;   //position of each unit in `ord`; NA for units not in it
  Rcpp::IntegerVector first; //first position in `ord` of each stratum
  Rcpp::IntegerVector last;  //last position in `ord` of each stratum
  int code_min = 0;          //stratum code of `first[0]` and `last[0]`
  int na_stratum = 0;        //index in `first` and `last` of the units with a missing code

  //Index in `first` and `last` of the stratum with code `code`, or -1 when no unit in
  //`ord` could have that code
  int stratum(int code) const {
    if (code == NA_INTEGER) {
      return na_stratum;
    }

    if (code < code_min || code - code_min >= na_stratum) {
      return -1;
    }

    return code - code_min;
  }
};

//`exact` holds the stratum code of every unit; `units` lists the units to arrange,
//in the order each stratum's units are to keep
ExactOrder make_exact_order(const Rcpp::IntegerVector& exact,
                            const Rcpp::IntegerVector& units);

//How find_control_mahcovs() and find_control_mat() treat `exact` strata. When `use`
//is set, only the treated unit's stratum is searched, in `order`. When `local` is
//also set, the search runs as a separate match of that stratum would run it rather
//than as a match of the whole sample would; the two give the same matches except in
//how they break ties between equally close controls. See nn_matchC_mahcovs().
struct StrataScan {
  bool use = false;
  bool local = false;
  Rcpp::IntegerVector strata; //stratum code of each unit
  ExactOrder order;
};

//The column of `mah_covs` to sort the units on for find_control_mahcovs() and the
//caliper on that column's scale: the first column that is an affine transformation
//of a caliper covariate with a nonnegative caliper, judged on the units in `rows`,
//or on all units when `rows` is empty. Column 0 and no caliper when there is none.
std::pair<int, double> find_match_var(const Rcpp::NumericMatrix& mah_covs,
                                      const Rcpp::NumericMatrix& caliper_covs_mat,
                                      const Rcpp::NumericVector& caliper_covs,
                                      const Rcpp::IntegerVector& rows);

//Sets up `scan`, the matching variable `match_var`, and the order the units are
//searched in (`ind_d_ord` and its inverse `match_d_ord`) for find_control_mahcovs().
//Returns the caliper on the matching variable's scale: one value for the whole
//sample, or with `local` one for each stratum, indexed by `scan.order.stratum()`.
std::vector<double> set_up_mahcovs_scan(StrataScan& scan,
                                        Rcpp::NumericVector& match_var,
                                        Rcpp::IntegerVector& ind_d_ord,
                                        Rcpp::IntegerVector& match_d_ord,
                                        const Rcpp::NumericMatrix& mah_covs,
                                        const Rcpp::NumericMatrix& caliper_covs_mat,
                                        const Rcpp::NumericVector& caliper_covs,
                                        const Rcpp::Nullable<Rcpp::IntegerVector>& strata_,
                                        bool local,
                                        Rcpp::Function& o);

//The caliper on the matching variable's scale for treated unit `t_id`, from the
//values set_up_mahcovs_scan() returns
inline double caliper_on_match_var(const std::vector<double>& match_var_caliper,
                                   const StrataScan& scan,
                                   int t_id) {
  return scan.local ? match_var_caliper[scan.order.stratum(scan.strata[t_id])] : match_var_caliper[0];
}

std::vector<int> find_control_vec(int t_id,
                                  const Rcpp::IntegerVector& ind_d_ord,
                                  const Rcpp::IntegerVector& match_d_ord,
                                  const Rcpp::IntegerVector& treat,
                                  const Rcpp::NumericVector& distance,
                                  const Rcpp::LogicalVector& eligible,
                                  int gi,
                                  int r,
                                  const Rcpp::IntegerVector& mm_rowi_,
                                  int ncc,
                                  const Rcpp::NumericMatrix& caliper_covs_mat,
                                  const Rcpp::NumericVector& caliper_covs,
                                  double caliper_dist,
                                  bool use_exact,
                                  const Rcpp::IntegerVector& exact,
                                  int aenc,
                                  const Rcpp::IntegerMatrix& antiexact_covs,
                                  const Rcpp::IntegerVector& first_control,
                                  const Rcpp::IntegerVector& last_control,
                                  const ExactOrder& exact_order,
                                  int ratio = 1,
                                  int prev_start = -1);

std::vector<int> find_control_mahcovs(int t_id,
                                      const Rcpp::IntegerVector& ind_d_ord,
                                      const Rcpp::IntegerVector& match_d_ord,
                                      const Rcpp::NumericVector& match_var,
                                      double match_var_caliper,
                                      const Rcpp::IntegerVector& treat,
                                      const Rcpp::NumericVector& distance,
                                      const Rcpp::LogicalVector& eligible,
                                      int gi,
                                      int r,
                                      const Rcpp::IntegerVector& mm_rowi,
                                      const Rcpp::NumericMatrix& mah_covs,
                                      int ncc,
                                      const Rcpp::NumericMatrix& caliper_covs_mat,
                                      const Rcpp::NumericVector& caliper_covs,
                                      bool use_caliper_dist,
                                      double caliper_dist,
                                      bool use_exact,
                                      const Rcpp::IntegerVector& exact,
                                      int aenc,
                                      const Rcpp::IntegerMatrix& antiexact_covs,
                                      const StrataScan& scan,
                                      int ratio = 1);

//Row `t_row` of `distance_mat` holds the treated unit's distances to the controls;
//`ind_match` gives each control's column
std::vector<int> find_control_mat(int t_id,
                                  const Rcpp::IntegerVector& treat,
                                  const Rcpp::IntegerVector& ind_non_focal,
                                  const Rcpp::IntegerVector& ind_match,
                                  const Rcpp::NumericMatrix& distance_mat,
                                  int t_row,
                                  const Rcpp::LogicalVector& eligible,
                                  int gi,
                                  int r,
                                  const Rcpp::IntegerVector& mm_rowi,
                                  int ncc,
                                  const Rcpp::NumericMatrix& caliper_covs_mat,
                                  const Rcpp::NumericVector& caliper_covs,
                                  double caliper_dist,
                                  bool use_exact,
                                  const Rcpp::IntegerVector& exact,
                                  int aenc,
                                  const Rcpp::IntegerMatrix& antiexact_covs,
                                  const StrataScan& scan,
                                  int ratio = 1);

double euc_dist_sq(const Rcpp::NumericMatrix& x,
                   int i,
                   int j);

//`ids` is taken by value so the early returns can move it rather than copy
std::vector<int> take_closest(std::vector<int> ids,
                              const std::vector<double>& dists,
                              int ratio);

bool antiexact_okay(int aenc,
                    int i,
                    int j,
                    const Rcpp::IntegerMatrix& antiexact_covs);

bool caliper_covs_okay(int ncc,
                       int i,
                       int j,
                       const Rcpp::NumericMatrix& caliper_covs_mat,
                       const Rcpp::NumericVector& caliper_covs);

bool caliper_dist_okay(bool use_caliper_dist,
                       int i,
                       int j,
                       const Rcpp::NumericVector& distance,
                       double caliper_dist);

bool mm_okay(int r,
             int i,
             const Rcpp::IntegerVector& mm_rowi);

bool exact_okay(bool use_exact,
                int i,
                int j,
                const Rcpp::IntegerVector& exact);

double max_finite(const Rcpp::NumericVector& x);

double min_finite(const Rcpp::NumericVector& x);

//`first_control` and `last_control` are taken by value on purpose: the copies share
//the caller's SEXP, which is how the updates below reach the caller.
void update_first_and_last_control(Rcpp::IntegerVector first_control,
                                   Rcpp::IntegerVector last_control,
                                   const Rcpp::IntegerVector& ind_d_ord,
                                   const Rcpp::LogicalVector& eligible,
                                   const Rcpp::IntegerVector& treat,
                                   int gi);

double get_affine_transformation(const Rcpp::NumericVector& x,
                                 const Rcpp::NumericVector& y,
                                 double tol = 1e-9);

#endif
