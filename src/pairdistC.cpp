#include "internal.h"
using namespace Rcpp;

//Sum of |x_i - x_j| over the pairs of a treated and a control unit among the units
//in `units`, which are sorted here on `x`. Each unit differs from every unit of the
//other group before it by its own value minus theirs, so the total is the sum over
//units of (number of other-group units so far) * value - (their sum so far): one
//pass rather than a visit to every pair. Values are taken relative to the smallest
//so that a large common offset in `x` does not cost precision.
static long double sorted_pair_sum(std::vector<std::pair<double, bool>>& units) {
  //Sorting is undefined with a `NaN` in the range, and any pair involving one is
  //missing anyway, as it is when the pairs are summed directly
  for (const auto& u : units) {
    if (std::isnan(u.first)) {
      return NA_REAL;
    }
  }

  std::sort(units.begin(), units.end());

  double x0 = units[0].first;

  long double sum1 = 0, sum0 = 0, total = 0;
  double n1 = 0, n0 = 0;

  for (const auto& u : units) {
    double v = u.first - x0;

    if (u.second) {
      total += n0 * v - sum0;
      n1++;
      sum1 += v;
    }
    else {
      total += n1 * v - sum1;
      n0++;
      sum0 += v;
    }
  }

  return total;
}

//Mean absolute difference in `x` over all pairs of a treated and a control unit in the
//same subclass. Units with a missing subclass are ignored.
//
//A subclass with few treated-control pairs, as after pair matching, is summed pair by
//pair. A larger one is sorted and summed with sorted_pair_sum(), whose time grows
//with the subclass's size rather than with its number of pairs, which is the square
//of its size.
// [[Rcpp::export]]
double pairdistsubC(const NumericVector& x,
                    const IntegerVector& t,
                    const IntegerVector& s) {

  //`base::order()`'s radix sort beats every C++ alternative measured here by 3-8x at
  //these sizes; see _dev/cpp-cleanup-notes.md. Looked up in the base environment
  //because `Function("order")` searches from the global environment, where a user
  //object of that name would mask it.
  Function ord = Environment::base_env()["order"];
  IntegerVector o = ord(s);
  o = o - 1;

  R_xlen_t n = sum(!is_na(s));

  const double max_direct_pairs = 256;

  long double total = 0;
  double n_pairs = 0;

  std::vector<std::pair<double, bool>> units;

  R_xlen_t i = 0;

  while (i < n) {
    //The subclass is o[i], ..., o[end - 1]
    R_xlen_t end = i + 1;
    while (end < n && s[o[end]] == s[o[i]]) {
      end++;
    }

    double n1 = 0;
    for (R_xlen_t a = i; a < end; a++) {
      if (t[o[a]] == 1) {
        n1++;
      }
    }

    double n0 = (end - i) - n1;

    if (n1 > 0 && n0 > 0) {
      if (n1 * n0 <= max_direct_pairs) {
        for (R_xlen_t a = i; a < end; a++) {
          if (t[o[a]] != 1) {
            continue;
          }

          for (R_xlen_t b = i; b < end; b++) {
            if (t[o[b]] == 1) {
              continue;
            }

            total += std::abs(x[o[a]] - x[o[b]]);
          }
        }
      }
      else {
        units.clear();

        for (R_xlen_t a = i; a < end; a++) {
          units.emplace_back(x[o[a]], t[o[a]] == 1);
        }

        total += sorted_pair_sum(units);
      }

      n_pairs += n1 * n0;
    }

    i = end;
  }

  if (n_pairs == 0) {
    return 0;
  }

  return static_cast<double>(total / n_pairs);
}
