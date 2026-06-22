# Snapshot tests for method = "nearest".
# These pin exact match.matrix results to detect changes in C++ matching algorithms.
# Complement test-method_nearest.R which tests structural correctness.

data("lalonde", package = "MatchIt")

# Fixed subset for distance matrix tests (50 treated + 50 control)
lalonde_sub <- lalonde[c(1:50, 186:235), ]

# ===== Baseline tests: one per C++ code path =====

test_that("baseline: PS vector, m.order='largest' (default)", {
  set.seed(12345)
  m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest")
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_snapshot_value(m$match.matrix, style = "json2")
})

test_that("baseline: Mahalanobis, m.order='data'", {
  set.seed(12345)
  m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest",
               mahvars = ~ age + educ + re74 + re75,
               m.order = "data")
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_snapshot_value(m$match.matrix, style = "json2")
})

test_that("baseline: distance matrix, m.order='data'", {
  set.seed(12345)
  d <- scaled_euclidean_dist(treat ~ age + educ + re74 + re75, data = lalonde_sub)
  m <- matchit(treat ~ age + educ + re74 + re75, data = lalonde_sub,
               distance = d, m.order = "data")
  expect_good_matchit(m, expect_distance = FALSE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_snapshot_value(m$match.matrix, style = "json2")
})

test_that("baseline: PS vector, m.order='closest'", {
  set.seed(12345)
  m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest",
               m.order = "closest")
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_snapshot_value(m$match.matrix, style = "json2")
})

test_that("baseline: Mahalanobis, m.order='closest'", {
  set.seed(12345)
  m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest",
               mahvars = ~ age + educ + re74 + re75,
               m.order = "closest")
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_snapshot_value(m$match.matrix, style = "json2")
})

test_that("baseline: distance matrix, m.order='closest'", {
  set.seed(12345)
  d <- scaled_euclidean_dist(treat ~ age + educ + re74 + re75, data = lalonde_sub)
  m <- matchit(treat ~ age + educ + re74 + re75, data = lalonde_sub,
               distance = d, m.order = "closest")
  expect_good_matchit(m, expect_distance = FALSE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_snapshot_value(m$match.matrix, style = "json2")
})

# ===== m.order variants =====

test_that("PS vector, m.order='data'", {
  set.seed(12345)
  m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest",
               m.order = "data")
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_snapshot_value(m$match.matrix, style = "json2")
})

test_that("PS vector, m.order='random'", {
  set.seed(12345)
  m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest",
               m.order = "random")
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_snapshot_value(m$match.matrix, style = "json2")
})

test_that("PS vector, m.order='farthest' (close=FALSE)", {
  set.seed(12345)
  m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest",
               m.order = "farthest")
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_snapshot_value(m$match.matrix, style = "json2")
})

test_that("Mahalanobis, m.order='random' (mahcovs path)", {
  set.seed(12345)
  m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest",
               mahvars = ~ age + educ + re74 + re75,
               m.order = "random")
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_snapshot_value(m$match.matrix, style = "json2")
})

test_that("Mahalanobis, m.order='farthest' (mahcovs_closest close=FALSE)", {
  set.seed(12345)
  m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest",
               mahvars = ~ age + educ + re74 + re75,
               m.order = "farthest")
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_snapshot_value(m$match.matrix, style = "json2")
})

test_that("distance matrix, m.order='random'", {
  set.seed(12345)
  d <- scaled_euclidean_dist(treat ~ age + educ + re74 + re75, data = lalonde_sub)
  m <- matchit(treat ~ age + educ + re74 + re75, data = lalonde_sub,
               distance = d, m.order = "random")
  expect_good_matchit(m, expect_distance = FALSE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_snapshot_value(m$match.matrix, style = "json2")
})

test_that("distance matrix, m.order='farthest' (distmat_closest close=FALSE)", {
  set.seed(12345)
  d <- scaled_euclidean_dist(treat ~ age + educ + re74 + re75, data = lalonde_sub)
  m <- matchit(treat ~ age + educ + re74 + re75, data = lalonde_sub,
               distance = d, m.order = "farthest")
  expect_good_matchit(m, expect_distance = FALSE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_snapshot_value(m$match.matrix, style = "json2")
})

# ===== ratio + replacement =====

test_that("ratio=3, replace=FALSE (pool depletion)", {
  set.seed(12345)
  expect_warning(
    m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
                 data = lalonde, method = "nearest",
                 ratio = 3, replace = FALSE),
    "Not all treated units will get 3 matches"
  )
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 3L)
  expect_snapshot_value(m$match.matrix, style = "json2")
})

test_that("ratio=3, replace=TRUE (no pool depletion)", {
  set.seed(12345)
  m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest",
               ratio = 3, replace = TRUE)
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = FALSE, ratio = 3L, replace = TRUE)
  expect_snapshot_value(m$match.matrix, style = "json2")
})

# ===== ratio + caliper =====

test_that("ratio=2, positive caliper (pool restriction)", {
  set.seed(12345)
  m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest",
               ratio = 2, caliper = 0.1, std.caliper = FALSE)
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 2L)
  expect_snapshot_value(m$match.matrix, style = "json2")
})

test_that("ratio=2, negative caliper (anti-caliper)", {
  set.seed(12345)
  m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest",
               ratio = 2, caliper = -0.05, std.caliper = FALSE)
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 2L)
  expect_snapshot_value(m$match.matrix, style = "json2")
})

# ===== caliper + replacement =====

test_that("caliper + replace=TRUE + ratio=2 (reuse within caliper)", {
  set.seed(12345)
  m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest",
               caliper = 0.1, std.caliper = FALSE,
               replace = TRUE, ratio = 2)
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = FALSE, ratio = 2L, replace = TRUE)
  expect_snapshot_value(m$match.matrix, style = "json2")
})

# ===== m.order + caliper + no replacement =====

test_that("m.order='largest' + caliper + replace=FALSE", {
  set.seed(12345)
  m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest",
               m.order = "largest", caliper = 0.1, std.caliper = FALSE,
               replace = FALSE)
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_snapshot_value(m$match.matrix, style = "json2")
})

test_that("m.order='smallest' + caliper + replace=FALSE", {
  set.seed(12345)
  m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest",
               m.order = "smallest", caliper = 0.1, std.caliper = FALSE,
               replace = FALSE)
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_snapshot_value(m$match.matrix, style = "json2")
})

# ===== exact + other constraints =====

test_that("exact + ratio=2 (within-stratum depletion)", {
  set.seed(12345)
  expect_warning(
    expect_warning(
      m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
                   data = lalonde, method = "nearest",
                   exact = ~ race, ratio = 2),
      "Fewer control units than treated units"
    ),
    "Not all treated units will get 2 matches"
  )
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 2L)
  expect_snapshot_value(m$match.matrix, style = "json2")
})

test_that("exact + caliper (double constraint)", {
  set.seed(12345)
  expect_warning(
    m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
                 data = lalonde, method = "nearest",
                 exact = ~ race, caliper = 0.2, std.caliper = FALSE),
    "Fewer control units than treated units"
  )
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_snapshot_value(m$match.matrix, style = "json2")
})

test_that("exact + antiexact (inclusion + exclusion)", {
  set.seed(12345)
  expect_warning(
    m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
                 data = lalonde, method = "nearest",
                 exact = ~ race, antiexact = ~ married),
    "Fewer control units than treated units"
  )
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_snapshot_value(m$match.matrix, style = "json2")
})

test_that("exact + replace=TRUE + ratio=2 (reuse within strata)", {
  set.seed(12345)
  m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest",
               exact = ~ race, replace = TRUE, ratio = 2)
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = FALSE, ratio = 2L, replace = TRUE)
  expect_snapshot_value(m$match.matrix, style = "json2")
})

# ===== mahvars + caliper =====

test_that("mahvars + distance caliper (Mahalanobis match with PS caliper)", {
  set.seed(12345)
  m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest",
               mahvars = ~ age + educ + re74 + re75,
               caliper = 0.2, std.caliper = FALSE)
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_snapshot_value(m$match.matrix, style = "json2")
})

test_that("mahvars + exact + m.order='closest' (three-way)", {
  set.seed(12345)
  expect_warning(
    m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
                 data = lalonde, method = "nearest",
                 mahvars = ~ age + educ + re74 + re75,
                 exact = ~ race, m.order = "closest"),
    "Fewer control units than treated units"
  )
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_snapshot_value(m$match.matrix, style = "json2")
})

# ===== mahvars + replace =====

test_that("mahvars + replace=TRUE (mahcovs path, reuse allowed)", {
  set.seed(12345)
  m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest",
               mahvars = ~ age + educ + re74 + re75,
               replace = TRUE)
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = FALSE, ratio = 1L, replace = TRUE)
  expect_snapshot_value(m$match.matrix, style = "json2")
})

test_that("mahvars + replace=TRUE + ratio=2 (mahcovs, reuse, multi-match)", {
  set.seed(12345)
  m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest",
               mahvars = ~ age + educ + re74 + re75,
               replace = TRUE, ratio = 2)
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = FALSE, ratio = 2L, replace = TRUE)
  expect_snapshot_value(m$match.matrix, style = "json2")
})

# ===== mahvars + antiexact =====

test_that("mahvars + antiexact (mahcovs path with antiexact)", {
  set.seed(12345)
  m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest",
               mahvars = ~ age + educ + re74 + re75,
               antiexact = ~ married)
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_snapshot_value(m$match.matrix, style = "json2")
})

# ===== reuse.max =====

test_that("reuse.max=3 + ratio=2 (bounded replacement)", {
  set.seed(12345)
  m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest",
               ratio = 2, reuse.max = 3)
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = FALSE, ratio = 2L,
                      replace = structure(TRUE, reuse.max = 3))
  expect_snapshot_value(m$match.matrix, style = "json2")
})

test_that("reuse.max=2 + ratio=2 + caliper (bounded replacement + caliper)", {
  set.seed(12345)
  m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest",
               ratio = 2, reuse.max = 2,
               caliper = 0.2, std.caliper = FALSE)
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = FALSE, ratio = 2L,
                      replace = structure(TRUE, reuse.max = 2))
  expect_snapshot_value(m$match.matrix, style = "json2")
})

test_that("mahvars + reuse.max=3 + ratio=2 (mahcovs bounded replacement)", {
  set.seed(12345)
  m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest",
               mahvars = ~ age + educ + re74 + re75,
               ratio = 2, reuse.max = 3)
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = FALSE, ratio = 2L,
                      replace = structure(TRUE, reuse.max = 3))
  expect_snapshot_value(m$match.matrix, style = "json2")
})

test_that("mahvars + reuse.max=2 + m.order='closest' (mahcovs_closest bounded)", {
  set.seed(12345)
  m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest",
               mahvars = ~ age + educ + re74 + re75,
               ratio = 2, reuse.max = 2, m.order = "closest")
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = FALSE, ratio = 2L,
                      replace = structure(TRUE, reuse.max = 2))
  expect_snapshot_value(m$match.matrix, style = "json2")
})

test_that("distmat + reuse.max=3 + ratio=2 (distmat bounded replacement)", {
  set.seed(12345)
  d <- scaled_euclidean_dist(treat ~ age + educ + re74 + re75, data = lalonde_sub)
  m <- matchit(treat ~ age + educ + re74 + re75, data = lalonde_sub,
               distance = d, ratio = 2, reuse.max = 3)
  expect_good_matchit(m, expect_distance = FALSE, expect_match.matrix = TRUE,
                      expect_subclass = FALSE, ratio = 2L,
                      replace = structure(TRUE, reuse.max = 3))
  expect_snapshot_value(m$match.matrix, style = "json2")
})

# ===== variable ratio =====

test_that("variable ratio + caliper (min/max with restriction)", {
  set.seed(12345)
  m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest",
               ratio = 2, min.controls = 1, max.controls = 4,
               caliper = 0.2, std.caliper = FALSE)
  ratio_attr <- structure(2L, min.controls = 1, max.controls = 4)
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = ratio_attr)
  expect_snapshot_value(m$match.matrix, style = "json2")
})

test_that("variable ratio + exact (within strata)", {
  set.seed(12345)
  expect_warning(
    expect_warning(
      m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
                   data = lalonde, method = "nearest",
                   ratio = 2, min.controls = 1, max.controls = 4,
                   exact = ~ race),
      "Fewer control units than treated units"
    ),
    "Not enough control units"
  )
  ratio_attr <- structure(2L, min.controls = 1, max.controls = 4)
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = ratio_attr)
  expect_snapshot_value(m$match.matrix, style = "json2")
})

test_that("variable ratio baseline (min/max.controls, PS vector)", {
  set.seed(12345)
  m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest",
               ratio = 2, min.controls = 1, max.controls = 4)
  ratio_attr <- structure(2L, min.controls = 1, max.controls = 4)
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = ratio_attr)
  expect_snapshot_value(m$match.matrix, style = "json2")
})

# ===== estimand = "ATC" =====

test_that("estimand='ATC' baseline (PS vector, flipped focal)", {
  set.seed(12345)
  expect_warning(
    m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
                 data = lalonde, method = "nearest",
                 estimand = "ATC"),
    "Fewer treated units than control units"
  )
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_snapshot_value(m$match.matrix, style = "json2")
})

test_that("estimand='ATC' + mahvars (mahcovs path, ATC focal)", {
  set.seed(12345)
  expect_warning(
    m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
                 data = lalonde, method = "nearest",
                 mahvars = ~ age + educ + re74 + re75,
                 estimand = "ATC"),
    "Fewer treated units than control units"
  )
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_snapshot_value(m$match.matrix, style = "json2")
})

test_that("estimand='ATC' + ratio=2 + exact (flipped focal)", {
  set.seed(12345)
  expect_warning(
    expect_warning(
      m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
                   data = lalonde, method = "nearest",
                   estimand = "ATC", ratio = 2, exact = ~ race),
      "Fewer treated units than control units"
    ),
    "Not all control units will get 2 matches"
  )
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 2L)
  expect_snapshot_value(m$match.matrix, style = "json2")
})

# ===== multiple calipers =====

test_that("covariate caliper + distance caliper (both simultaneously)", {
  set.seed(12345)
  m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest",
               caliper = c(.1, age = 2), std.caliper = FALSE)
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_snapshot_value(m$match.matrix, style = "json2")
})

test_that("covariate caliper + antiexact + m.order='closest'", {
  set.seed(12345)
  m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest",
               caliper = c(age = 5), std.caliper = FALSE,
               antiexact = ~ married, m.order = "closest")
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_snapshot_value(m$match.matrix, style = "json2")
})

# ===== std.caliper = TRUE =====

test_that("standardized caliper (std.caliper=TRUE, PS vector)", {
  set.seed(12345)
  m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest",
               caliper = 0.25, std.caliper = TRUE)
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_snapshot_value(m$match.matrix, style = "json2")
})

test_that("standardized covariate caliper (std.caliper=TRUE on age)", {
  set.seed(12345)
  m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest",
               caliper = c(age = 1), std.caliper = TRUE)
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_snapshot_value(m$match.matrix, style = "json2")
})

# ===== antiexact without exact =====

test_that("antiexact alone, PS vector", {
  set.seed(12345)
  m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest",
               antiexact = ~ married)
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_snapshot_value(m$match.matrix, style = "json2")
})

test_that("antiexact + m.order='closest' (vec_closest path)", {
  set.seed(12345)
  m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest",
               antiexact = ~ married, m.order = "closest")
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_snapshot_value(m$match.matrix, style = "json2")
})

# ===== discard =====

test_that("discard logical vector, PS vector path", {
  set.seed(12345)
  dis <- lalonde$re74 > 15000
  m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest",
               discard = dis)
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_snapshot_value(m$match.matrix, style = "json2")
})

test_that("discard + ratio=2 + m.order='closest' (vec_closest path with discards)", {
  set.seed(12345)
  dis <- lalonde$re74 > 10000
  m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest",
               discard = dis, ratio = 2, m.order = "closest")
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 2L, replace = FALSE)
  expect_snapshot_value(m$match.matrix, style = "json2")
})

test_that("discard + mahvars (mahcovs path with discards)", {
  set.seed(12345)
  dis <- lalonde$re74 > 15000
  m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest",
               mahvars = ~ age + educ + re74 + re75,
               discard = dis)
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_snapshot_value(m$match.matrix, style = "json2")
})

# ===== unit.id =====

test_that("unit.id with replacement=FALSE (clustered units, PS vector)", {
  set.seed(12345)
  m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest",
               unit.id = ~ age)
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_snapshot_value(m$match.matrix, style = "json2")
})

test_that("unit.id + m.order='closest' (vec_closest path with unit.id)", {
  set.seed(12345)
  expect_warning(
    m <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
                 data = lalonde, method = "nearest",
                 unit.id = ~ age, m.order = "closest", ratio = 2),
    "Not all treated units will get 2 matches"
  )
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 2L)
  expect_snapshot_value(m$match.matrix, style = "json2")
})

# ===== distance matrix + constraints =====

test_that("distmat + exact + ratio=2 (stratum loop + multi-match)", {
  set.seed(12345)
  d <- scaled_euclidean_dist(treat ~ age + educ + re74 + re75, data = lalonde_sub)
  expect_warning(
    expect_warning(
      m <- matchit(treat ~ age + educ + re74 + re75, data = lalonde_sub,
                   distance = d, exact = ~ race, ratio = 2),
      "Fewer control units than treated units"
    ),
    "Not all treated units will get 2 matches"
  )
  expect_good_matchit(m, expect_distance = FALSE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 2L)
  expect_snapshot_value(m$match.matrix, style = "json2")
})

test_that("distmat + m.order='closest' + replace=TRUE + ratio=2", {
  set.seed(12345)
  d <- scaled_euclidean_dist(treat ~ age + educ + re74 + re75, data = lalonde_sub)
  m <- matchit(treat ~ age + educ + re74 + re75, data = lalonde_sub,
               distance = d, m.order = "closest",
               replace = TRUE, ratio = 2)
  expect_good_matchit(m, expect_distance = FALSE, expect_match.matrix = TRUE,
                      expect_subclass = FALSE, ratio = 2L, replace = TRUE)
  expect_snapshot_value(m$match.matrix, style = "json2")
})

test_that("distmat + caliper (distmat path with caliper constraint)", {
  set.seed(12345)
  d <- scaled_euclidean_dist(treat ~ age + educ + re74 + re75, data = lalonde_sub)
  m <- matchit(treat ~ age + educ + re74 + re75, data = lalonde_sub,
               distance = d, caliper = c(age = 3), std.caliper = FALSE)
  expect_good_matchit(m, expect_distance = FALSE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_snapshot_value(m$match.matrix, style = "json2")
})

test_that("distmat + antiexact (distmat path with antiexact constraint)", {
  set.seed(12345)
  d <- scaled_euclidean_dist(treat ~ age + educ + re74 + re75, data = lalonde_sub)
  m <- matchit(treat ~ age + educ + re74 + re75, data = lalonde_sub,
               distance = d, antiexact = ~ married)
  expect_good_matchit(m, expect_distance = FALSE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_snapshot_value(m$match.matrix, style = "json2")
})

# ===== full Mahalanobis distance =====

test_that("full Mahalanobis (distance='mahalanobis')", {
  set.seed(12345)
  m <- matchit(treat ~ age + educ + re74 + re75,
               data = lalonde, method = "nearest",
               distance = "mahalanobis")
  expect_good_matchit(m, expect_distance = FALSE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_snapshot_value(m$match.matrix, style = "json2")
})

test_that("full Mahalanobis + m.order='closest'", {
  set.seed(12345)
  m <- matchit(treat ~ age + educ + re74 + re75,
               data = lalonde, method = "nearest",
               distance = "mahalanobis", m.order = "closest")
  expect_good_matchit(m, expect_distance = FALSE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L)
  expect_snapshot_value(m$match.matrix, style = "json2")
})
