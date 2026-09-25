test_that("distance vector, mah vars, and distance matrix yield identical results", {
  set.seed(1234)
  n <- 1e3
  p <- runif(n, 0, .4)
  x <- runif(n)
  g <- sample(1:5, n, TRUE)
  a <- rbinom(n, 1, p)
  u <- 1:n; u[a == 0] <- sample(u[a == 0][1:round(sum(a == 0)/5)], sum(a == 0), replace = TRUE)
  dis <- as.logical(rbinom(n, 1, .1))
  d <- data.frame(p, x, a, g, u, dis)
  d$p_ <- d$p

  dd <- euclidean_dist(a ~ p, data = d)

  test_all <- function(..., which = 1:4) {

    M <- list()
    if (any(which == 1)) {
      matchit_try({
        m <- matchit(a ~ p + p_, data = d,
                     distance = d$p,
                     ...)
      }, dont_warn_if = "Fewer control")
      expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                          expect_subclass = !m$info$replace, replace = m$info$replace,
                          ratio = m$info$ratio)
      M <- c(M, list(m))
    }
    if (any(which == 2)) {
      matchit_try({
        m <- matchit(a ~ p + p_, data = d,
                     distance = "euclidean",
                     ...)
      }, dont_warn_if = "Fewer control")
      expect_good_matchit(m, expect_distance = FALSE, expect_match.matrix = TRUE,
                          expect_subclass = !m$info$replace, replace = m$info$replace,
                          ratio = m$info$ratio)
      M <- c(M, list(m))
    }
    if (any(which == 3)) {
      matchit_try({
        m <- matchit(a ~ p + p_, data = d,
                     distance = d$p,
                     mahvars = ~p + p_,
                     ...)
      }, dont_warn_if = "Fewer control")
      expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                          expect_subclass = !m$info$replace, replace = m$info$replace,
                          ratio = m$info$ratio)
      M <- c(M, list(m))
    }
    if (any(which == 4)) {
      matchit_try({
        m <- matchit(a ~ p + p_, data = d,
                     distance = dd,
                     ...)
      }, dont_warn_if = "Fewer control")
      expect_good_matchit(m, expect_distance = FALSE, expect_match.matrix = TRUE,
                          expect_subclass = !m$info$replace, replace = m$info$replace,
                          ratio = m$info$ratio)
      M <- c(M, list(m))
    }

    all(unlist(lapply(M[-1], function(m) isTRUE(all.equal(M[[1]]$match.matrix,
                                                          m$match.matrix)))))
  }

  expect_true(test_all(m.order = "data"))
  expect_true(test_all(m.order = "closest"))
  expect_true(test_all(m.order = "farthest"))
  expect_true(test_all(m.order = "largest", which = c(1, 3)))
  expect_true(test_all(m.order = "smallest", which = c(1, 3)))

  expect_true(test_all(m.order = "data", ratio = 2))
  expect_true(test_all(m.order = "closest", ratio = 2))

  expect_true(test_all(m.order = "data", ratio = 2, max.controls = 3, which = c(1, 3)))
  expect_true(test_all(m.order = "closest", ratio = 2, max.controls = 3, which = c(1, 3)))
  expect_true(test_all(m.order = "largest", ratio = 2, max.controls = 3, which = c(1, 3)))

  expect_true(test_all(m.order = "data", ratio = 2, replace = TRUE))
  expect_true(test_all(m.order = "closest", ratio = 2, replace = TRUE))

  expect_true(test_all(m.order = "data", ratio = 2, reuse.max = 3))
  expect_true(test_all(m.order = "closest", ratio = 2, reuse.max = 3))

  expect_true(test_all(m.order = "data", ratio = 2, caliper = .001, std.caliper = FALSE, which = c(1, 3)))
  expect_true(test_all(m.order = "closest", ratio = 2, caliper = .001, std.caliper = FALSE, which = c(1, 3)))
  expect_true(test_all(m.order = "largest", ratio = 2, caliper = .001, std.caliper = FALSE, which = c(1, 3)))

  expect_true(test_all(m.order = "data", ratio = 2, caliper = -.001, std.caliper = FALSE, which = c(1, 3)))
  expect_true(test_all(m.order = "closest", ratio = 2, caliper = -.001, std.caliper = FALSE, which = c(1, 3)))
  expect_true(test_all(m.order = "largest", ratio = 2, caliper = -.001, std.caliper = FALSE, which = c(1, 3)))

  expect_true(test_all(m.order = "data", ratio = 2, caliper = c(p = .001), std.caliper = FALSE))
  expect_true(test_all(m.order = "closest", ratio = 2, caliper = c(p = .001), std.caliper = FALSE))

  expect_true(test_all(m.order = "data", ratio = 2, caliper = c(p = -.001), std.caliper = FALSE))
  expect_true(test_all(m.order = "closest", ratio = 2, caliper = c(p = -.001), std.caliper = FALSE))

  expect_true(test_all(m.order = "data", ratio = 2, caliper = c(p = .001), std.caliper = FALSE, reuse.max = 3))
  expect_true(test_all(m.order = "closest", ratio = 2, caliper = c(p = .001), std.caliper = FALSE, reuse.max = 3))

  expect_true(test_all(m.order = "data", ratio = 2, caliper = c(p = -.001), std.caliper = FALSE, reuse.max = 3))
  expect_true(test_all(m.order = "closest", ratio = 2, caliper = c(p = -.001), std.caliper = FALSE, reuse.max = 3))

  expect_true(test_all(m.order = "data", ratio = 2, exact = ~g))
  expect_true(test_all(m.order = "closest", ratio = 2, exact = ~g))

  expect_true(test_all(m.order = "data", ratio = 2, exact = ~g, replace = TRUE))
  expect_true(test_all(m.order = "closest", ratio = 2, exact = ~g, replace = TRUE))

  expect_true(test_all(m.order = "data", ratio = 2, antiexact = ~g))
  expect_true(test_all(m.order = "closest", ratio = 2, antiexact = ~g))

  expect_true(test_all(m.order = "data", ratio = 2, antiexact = ~g, replace = TRUE))
  expect_true(test_all(m.order = "closest", ratio = 2, antiexact = ~g, replace = TRUE))

  expect_true(test_all(m.order = "data", ratio = 2, discard = dis))
  expect_true(test_all(m.order = "closest", ratio = 2, discard = dis))

  expect_true(test_all(m.order = "data", ratio = 2, unit.id = ~u))
  expect_true(test_all(m.order = "closest", ratio = 2, unit.id = ~u))

  expect_true(test_all(m.order = "data", ratio = 2, unit.id = ~u, reuse.max = 3))
  expect_true(test_all(m.order = "closest", ratio = 2, unit.id = ~u, reuse.max = 3))

  expect_true(test_all(m.order = "data", ratio = 2, unit.id = ~u, replace = TRUE))
  expect_true(test_all(m.order = "closest", ratio = 2, unit.id = ~u, replace = TRUE))
})

test_that("calipers work, positive and negative", {
  set.seed(1234)
  n <- 1e3
  p <- runif(n, 0, .4)
  x <- runif(n)
  g <- sample(1:5, n, TRUE)
  a <- rbinom(n, 1, p)
  u <- 1:n; u[a == 0] <- sample(u[a == 0][1:round(sum(a == 0)/5)], sum(a == 0), replace = TRUE)
  dis <- as.logical(rbinom(n, 1, .1))
  d <- data.frame(p, x, a, g, u, dis)

  #Positive calipers
  m <- matchit(a ~ p + x, data = d, distance = d$p, caliper = .001,
               std.caliper = FALSE, replace = TRUE, ratio = 3)

  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = !m$info$replace,
                      ratio = m$info$ratio)

  m <- matchit(a ~ p + x, data = d, distance = d$p, caliper = c(x = .001),
               std.caliper = FALSE, replace = TRUE, ratio = 3)

  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = !m$info$replace,
                      ratio = m$info$ratio)

  m <- matchit(a ~ p + x, data = d, distance = d$p, caliper = c(.02, x = .01),
               std.caliper = FALSE, replace = TRUE, ratio = 3)

  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = !m$info$replace,
                      ratio = m$info$ratio)

  #Negative calipers
  m <- matchit(a ~ p + x, data = d, distance = d$p, caliper = -.001,
               std.caliper = FALSE, replace = TRUE, ratio = 3)

  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = !m$info$replace,
                      ratio = m$info$ratio)

  m <- matchit(a ~ p + x, data = d, distance = d$p, caliper = c(x = -.001),
               std.caliper = FALSE, replace = TRUE, ratio = 3)

  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = !m$info$replace,
                      ratio = m$info$ratio)

  m <- matchit(a ~ p + x, data = d, distance = d$p, caliper = c(-.02, x = -.01),
               std.caliper = FALSE, replace = TRUE, ratio = 3)

  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = !m$info$replace,
                      ratio = m$info$ratio)

  m <- matchit(a ~ p + x, data = d, distance = d$p, caliper = c(-.02, x = .01),
               std.caliper = FALSE, replace = TRUE, ratio = 3)

  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = !m$info$replace,
                      ratio = m$info$ratio)

  m <- matchit(a ~ p + x, data = d, distance = d$p, caliper = c(.02, x = -.01),
               std.caliper = FALSE, replace = TRUE, ratio = 3)

  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = !m$info$replace,
                      ratio = m$info$ratio)
})

test_that("covariate calipers are respected in every `exact` stratum", {
  set.seed(1234)
  n <- 600
  g <- sample(LETTERS[1:10], n, TRUE)
  x <- rnorm(n)
  y <- rnorm(n)
  a <- rbinom(n, 1, .25)
  d <- data.frame(a, x, y, g)

  #`x` is both a matching variable and the caliper variable, so the caliper is
  #converted to the scale of the matching variables internally. That conversion
  #must not change the caliper itself, which is reused for each `exact` stratum.
  for (m.order in c("data", "random", "closest", "farthest")) {
    m <- matchit(a ~ x + y, data = d, distance = "mahalanobis", exact = ~g,
                 caliper = c(x = .5), std.caliper = FALSE, m.order = m.order)

    expect_good_matchit(m, expect_distance = FALSE, expect_match.matrix = TRUE,
                        expect_subclass = TRUE, ratio = 1)

    mm <- m$match.matrix

    caliper.diff <- vapply(seq_len(nrow(mm)), function(i) {
      ctrl <- na.omit(mm[i, ])

      if (length(ctrl) == 0L) {
        return(0)
      }

      max(abs(d[rownames(mm)[i], "x"] - d[ctrl, "x"]))
    }, numeric(1L))

    expect_true(all(caliper.diff <= .5))
  }
})

test_that("covariate calipers keep controls whose difference equals the caliper", {
  #When the distance is an affine transformation of a caliper covariate, the search
  #for controls stops on the caliper restated on the distance's scale. Rounding in
  #that comparison used to stop it at controls whose difference on the covariate
  #equaled the caliper, or for a negative caliper barely exceeded it, which the
  #caliper itself accepts. With `m.order = "data"`, each treated unit must get the
  #closest control that is still available and within the caliper.
  closer_controls_skipped <- function(m, formula, data, caliper, exact = NULL) {
    X <- transform_covariates(formula, data = data, method = "scaled_euclidean",
                              treat = data$treat)
    v <- names(caliper)
    available <- rownames(data)[data$treat == 0]
    skipped <- character()

    for (i in rownames(m$match.matrix)) {
      diffs <- abs(data[i, v] - data[available, v])

      ok <- {
        if (caliper >= 0) diffs <= caliper
        else diffs > -caliper
      }

      if (is_not_null(exact)) {
        ok <- ok & data[available, exact] == data[i, exact]
      }

      within <- available[ok]
      dists <- sqrt(colSums((t(X[within, , drop = FALSE]) - X[i, ])^2))
      ctrl <- m$match.matrix[i, 1L]

      if (is.na(ctrl)) {
        if (is_not_null(within)) {
          skipped <- c(skipped, i)
        }

        next
      }

      if (sqrt(sum((X[ctrl, ] - X[i, ])^2)) > min(dists) + 1e-10) {
        skipped <- c(skipped, i)
      }

      available <- setdiff(available, ctrl)
    }

    skipped
  }

  #A Mahalanobis-type search, over the whole sample and within `exact` strata
  data("lalonde", package = "MatchIt", envir = environment())

  m <- matchit(treat ~ age + educ + married, data = lalonde,
               distance = "scaled_euclidean", caliper = c(educ = 1),
               std.caliper = FALSE, m.order = "data")

  expect_identical(unname(m$match.matrix["NSW4", 1L]), "PSID423")
  expect_identical(closer_controls_skipped(m, ~ age + educ + married, lalonde,
                                           c(educ = 1)),
                   character())

  m <- suppressWarnings({
    matchit(treat ~ age + educ + married, data = lalonde,
            distance = "scaled_euclidean", caliper = c(educ = 1), exact = ~ race,
            std.caliper = FALSE, m.order = "data")
  })

  expect_identical(closer_controls_skipped(m, ~ age + educ + married, lalonde,
                                           c(educ = 1), exact = "race"),
                   character())

  #A single covariate, matched as a distance vector, with a negative caliper
  set.seed(8)
  d <- data.frame(treat = rbinom(300, 1, .4), x = sample(0:20, 300, TRUE) / 10)
  rownames(d) <- paste0("u", seq_len(nrow(d)))

  m <- matchit(treat ~ x, data = d, distance = "scaled_euclidean",
               caliper = c(x = -.3), std.caliper = FALSE, m.order = "data")

  expect_identical(closer_controls_skipped(m, ~ x, d, c(x = -.3)),
                   character())
})

test_that("exact matching on many strata finds the same matches as each stratum alone", {
  #With `exact` and a distance vector, the search for a control scans only the
  #treated unit's stratum. Strata cannot share controls, so with no ties in the
  #distance, matching each stratum on its own must give the same matches.
  set.seed(4321)
  n <- 800
  d <- data.frame(p = runif(n), g = sample(1:40, n, TRUE))
  d$a <- rbinom(n, 1, .3)

  for (args in list(list(),
                    list(ratio = 2),
                    list(m.order = "closest"),
                    list(caliper = .01, std.caliper = FALSE),
                    list(replace = TRUE, ratio = 2))) {
    m <- suppressWarnings({
      do.call(matchit, c(list(a ~ p, data = d, distance = d$p, exact = ~g), args))
    })

    for (g in unique(d$g)) {
      d_g <- d[d$g == g, ]

      if (!all(0:1 %in% d_g$a)) {
        next
      }

      m_g <- tryCatch(suppressWarnings({
        do.call(matchit, c(list(a ~ p, data = d_g, distance = d_g$p), args))
      }), error = function(e) e)

      #A stratum where no treated unit has a control within the caliper
      if (inherits(m_g, "error")) {
        expect_match(conditionMessage(m_g), "No units were matched", fixed = TRUE)
        expect_true(all(is.na(m$match.matrix[rownames(d_g)[d_g$a == 1], ])))
        next
      }

      mm_g <- m_g$match.matrix
      expect_identical(m$match.matrix[rownames(mm_g), seq_len(ncol(mm_g)), drop = FALSE],
                       mm_g)
    }
  }
})

test_that("searching only the `exact` stratum finds the matches a search of the whole sample finds", {
  #With `unit.id` and a Mahalanobis distance or a distance matrix, the search for a
  #treated unit's controls visits only its stratum, but in the order a search of the
  #whole sample would reach them, which decides between equally close controls. The
  #covariates are discrete so that there are many such ties.
  set.seed(2468)
  n <- 400
  X <- matrix(sample(c(0, 1, 2, 3), 3 * n, TRUE), ncol = 3,
              dimnames = list(NULL, c("x1", "x2", "x3")))
  treat <- setNames(rbinom(n, 1, .35), paste0("u", seq_len(n)))
  n1 <- sum(treat == 1)
  g <- factor(sample(1:8, n, TRUE))
  uid <- factor(sample(1:300, n, TRUE))
  dm <- as.matrix(dist(X))[treat == 1, treat == 0]

  specs <- list(list(m.order = "data"),
                list(m.order = "random", ratio = 2L),
                list(m.order = "closest", ratio = 2L),
                list(m.order = "farthest"),
                list(m.order = "data", ratio = 3L, reuse.max = 2L),
                list(m.order = "data", unit.id = uid),
                list(m.order = "closest", ratio = 2L, unit.id = uid),
                list(m.order = "data", caliper.covs = c(x1 = 1)),
                list(m.order = "data", antiexactcovs = cbind(as.integer(X[, "x3"]))))

  for (s in specs) {
    for (d in c("mahcovs", "matrix")) {
      #`strata = NULL` searches the whole sample, checking `ex` for each pair
      match_in <- function(strata) {
        set.seed(1)
        nn_matchC_dispatch(treat, focal = 1L, ratio = rep.int(s$ratio %or% 1L, n1),
                           discarded = logical(n), reuse.max = s$reuse.max %or% 1L,
                           distance = NULL, distance_mat = if (d == "matrix") dm,
                           ex = g, caliper.dist = NULL, caliper.covs = s$caliper.covs,
                           caliper.covs.mat = if (is_not_null(s$caliper.covs)) {
                             X[, names(s$caliper.covs), drop = FALSE]
                           },
                           mahcovs = if (d == "mahcovs") X,
                           antiexactcovs = s$antiexactcovs, unit.id = s$unit.id,
                           m.order = s$m.order, verbose = FALSE, strata = strata)
      }

      expect_identical(match_in(g), match_in(NULL))
    }
  }
})

test_that("with unit.id, the exact-strata warning counts distinct control unit IDs", {
  #Stratum 1 has three treated units and four control units, but the controls have
  #only two unit IDs among them
  d <- data.frame(a = c(1, 1, 1, 0, 0, 0, 0, 1, 0, 0),
                  g = c(1, 1, 1, 1, 1, 1, 1, 2, 2, 2),
                  u = c(1, 2, 3, 4, 4, 5, 5, 6, 7, 8),
                  p = seq(.1, .9, length.out = 10))

  expect_no_unexpected_warning(matchit(a ~ p, data = d, distance = d$p, exact = ~g))

  expect_wrn(matchit(a ~ p, data = d, distance = d$p, exact = ~g, unit.id = ~u),
             "fewer control unit IDs than treated units in some `exact` strata")

  #A treated unit that shares an ID with a control does not remove that control's ID
  #from the count
  d2 <- data.frame(a = c(1, 1, 0, 0), g = 1, u = c(1, 2, 1, 2),
                   p = c(.2, .4, .3, .5))

  expect_no_unexpected_warning(matchit(a ~ p, data = d2, distance = d2$p, exact = ~g,
                                       unit.id = ~u))
})
