## Functions to calculate summary stats
bal1var <- function(xx, tt, ww = NULL, s.weights, subclass = NULL, mm = NULL,
                    s.d.denom = "treated", standardize = FALSE,
                    compute.pair.dist = TRUE) {

  un <- is_null(ww)
  bin.var <- all(xx == 0 | xx == 1)

  xsum <- rep.int(NA_real_, 7L)
  names(xsum) <- {
    if (standardize) c("Means Treated", "Means Control", "Std. Mean Diff.",
                       "Var. Ratio", "eCDF Mean", "eCDF Max", "Std. Pair Dist.")
    else c("Means Treated", "Means Control", "Mean Diff.",
           "Var. Ratio", "eQQ Mean", "eQQ Max", "Pair Dist.")
  }

  if (un) ww <- s.weights
  else ww <- ww * s.weights

  i1 <- which(tt == 1)
  i0 <- which(tt == 0)

  too.small <- sum(ww[i1] != 0) < 2L && sum(ww[i0] != 0) < 2L

  xsum["Means Treated"] <- wm(xx[i1], ww[i1], na.rm = TRUE)
  xsum["Means Control"] <- wm(xx[i0], ww[i0], na.rm = TRUE)

  mdiff <- xsum["Means Treated"] - xsum["Means Control"]

  if (standardize && abs(mdiff) > sqrt(.Machine$double.eps)) {
    if (!too.small) {
      std <- smd_denom(xx, tt, s.weights, s.d.denom, bin.var, i1 = i1, i0 = i0)

      xsum[3L] <- mdiff / std
      if (!un && compute.pair.dist) {
        xsum[7L] <- pair.dist(xx, tt, subclass, mm, std)
      }
    }
  }
  else {
    xsum[3L] <- mdiff

    if (!un && compute.pair.dist) {
      xsum[7L] <- pair.dist(xx, tt, subclass, mm)
    }
  }

  if (bin.var) {
    xsum[5L:6L] <- abs(mdiff)
  }
  else if (!too.small) {
    xsum["Var. Ratio"] <- wvar(xx[i1], bin.var, ww[i1]) / wvar(xx[i0], bin.var, ww[i0])

    qqmat <- qqsum(xx, tt, ww, standardize = standardize)
    xsum[5L:6L] <- qqmat[c("meandiff", "maxdiff")]
  }

  xsum
}

#Within-subclass balance statistics. `xx`, `tt`, and `s.weights` hold the units of one
#subclass only. `bin.var` and `std` describe the variable in the full sample: whether
#it is binary, and the denominator of its standardized mean difference from
#smd_denom() (SD from full sample, not within subclass). The caller computes them
#once per variable rather than once per subclass.
bal1var.subclass <- function(xx, tt, s.weights, bin.var, std = NULL,
                             standardize = FALSE) {

  xsum <- make_matrix(ncol = 6L, nrow = "Subclass")
  colnames(xsum) <- {
    if (standardize) c("Means Treated", "Means Control", "Std. Mean Diff.",
                       "Var. Ratio", "eCDF Mean", "eCDF Max")
    else c("Means Treated", "Means Control", "Mean Diff",
           "Var. Ratio", "eQQ Mean", "eQQ Max")
  }

  i1 <- which(tt == 1)
  i0 <- which(tt == 0)

  too.small <- length(i1) < 2L && length(i0) < 2L

  xsum["Subclass", "Means Treated"] <- wm(xx[i1], s.weights[i1], na.rm = TRUE)
  xsum["Subclass", "Means Control"] <- wm(xx[i0], s.weights[i0], na.rm = TRUE)

  mdiff <- xsum["Subclass", "Means Treated"] - xsum["Subclass", "Means Control"]

  if (standardize && abs(mdiff) > 1e-8) {
    if (!too.small) {
      xsum["Subclass", 3L] <- mdiff / std
    }
  }
  else {
    xsum["Subclass", 3L] <- mdiff
  }

  if (bin.var) {
    xsum["Subclass", 5L:6L] <- abs(mdiff)
  }
  else if (!too.small) {
    xsum["Subclass", "Var. Ratio"] <- wvar(xx[i1], bin.var, s.weights[i1]) / wvar(xx[i0], bin.var, s.weights[i0])

    qqall <- qqsum(xx, tt, s.weights, standardize = standardize)
    xsum["Subclass", 5L:6L] <- qqall[c("meandiff", "maxdiff")]
  }

  xsum
}

#Standardization factor for a standardized mean difference: `s.d.denom` itself when
#it is a number, and otherwise the standard deviation of `xx` in the group or groups
#it names, computed in the full sample with the sampling weights. The matching
#weights never enter it, so the unmatched, matched, and within-subclass statistics
#for a variable are all standardized by the same number. A caller that already has
#the indices of the treated and control units can pass them as `i1` and `i0`.
smd_denom <- function(xx, tt, s.weights, s.d.denom = "treated", bin.var = NULL,
                      i1 = which(tt == 1), i0 = which(tt == 0)) {
  if (is.numeric(s.d.denom)) {
    return(s.d.denom)
  }

  s.d.denom <- arg::match_arg(s.d.denom, c("treated", "control", "pooled"))

  std <- switch(s.d.denom,
                "treated" = sqrt(wvar(xx[i1], bin.var, s.weights[i1])),
                "control" = sqrt(wvar(xx[i0], bin.var, s.weights[i0])),
                "pooled" = pooled_sd(xx, tt, w = s.weights, bin.var = bin.var,
                                     contribution = "equal"))

  #Avoid divide by zero
  if (!is.finite(std) || std < sqrt(.Machine$double.eps)) {
    std <- pooled_sd(xx, tt, w = s.weights, bin.var = bin.var)
  }

  std
}

#Compute within-pair/subclass distances
pair.dist <- function(xx, tt, subclass = NULL, mm = NULL, std = NULL) {

  if (is_not_null(subclass)) {
    mpdiff <- pairdistsubC(as.numeric(xx), as.integer(tt),
                           as.integer(subclass))
  }
  else if (is_not_null(mm)) {
    names(xx) <- names(tt)
    xx_t <- xx[rownames(mm)]
    xx_c <- matrix(0, nrow = nrow(mm), ncol = ncol(mm))
    xx_c[] <- xx[mm]

    mpdiff <- mean(abs(xx_t - xx_c), na.rm = TRUE)
  }
  else {
    return(NA_real_)
  }

  if (is_not_null(std) && abs(mpdiff) > 1e-8) {
    return(mpdiff / std)
  }

  mpdiff
}

## Function for QQ summary stats
qqsum <- function(x, t, w = NULL, standardize = FALSE) {
  #x = variable, t = treat, w = weights

  #Names cannot affect the result, and copying them through every reordering and
  #subset below dominates the run time on large samples
  x <- unname(x)
  t <- unname(t)

  n.obs <- length(x)

  w <- {
    if (is_null(w)) rep.int(1, n.obs)
    else unname(w)
  }

  if (has_n_unique(x, 2L) && all(x == 0 | x == 1)) {
    t1 <- t == t[1L]
    #For binary variables, just difference in means. `wm()` divides by each group's
    #total weight, so unlike below the weights need no rescaling first.
    ediff <- abs(wm(x[t1], w[t1]) - wm(x[!t1], w[!t1]))

    return(c(meandiff = ediff, maxdiff = ediff))
  }

  w <- .make_sum_to_1(w, by = t)

  ord <- order(x)
  x_ord <- x[ord]
  w_ord <- w[ord]
  t_ord <- t[ord]

  t1 <- which(t_ord == t_ord[1L])

  if (standardize) {
    #Difference between ecdf of x for each group
    w_ord_ <- w_ord
    w_ord_[t1] <- -w_ord_[t1]
    ediff <- abs(cumsum(w_ord_))[c(diff1(x_ord) != 0, TRUE)]
  }
  else {
    #Horizontal distance of ecdf between groups
    #Need to interpolate larger group to be same size as smaller group

    u <- unique(x_ord)

    w1 <- w_ord[t1]
    w0 <- w_ord[-t1]

    x1 <- x_ord[t1][w1 > 0]
    x0 <- x_ord[-t1][w0 > 0]

    w1 <- w1[w1 > 0]
    w0 <- w0[w0 > 0]

    wn1 <- length(w1)
    wn0 <- length(w0)

    if (wn1 < wn0) {
      if (length(u) <= 5L) {
        x0probs <- vapply(u, function(u_) wm(x0 == u_, w0), numeric(1L))
        x0cumprobs <- c(0, .cumsum_prob(x0probs))
        x0 <- u[findInterval(.cumsum_prob(w1), x0cumprobs, rightmost.closed = TRUE)]
      }
      else {
        x0 <- approx(.cumsum_prob(w0), y = x0, xout = .cumsum_prob(w1), rule = 2,
                     method = "constant", ties = "ordered")$y
      }
    }
    else if (wn1 > wn0) {
      if (length(u) <= 5L) {
        x1probs <- vapply(u, function(u_) wm(x1 == u_, w1), numeric(1L))
        x1cumprobs <- c(0, .cumsum_prob(x1probs))
        x1 <- u[findInterval(.cumsum_prob(w0), x1cumprobs, rightmost.closed = TRUE)]
      }
      else {
        x1 <- approx(.cumsum_prob(w1), y = x1, xout = .cumsum_prob(w0), rule = 2,
                     method = "constant", ties = "ordered")$y
      }
    }

    ediff <- abs(x1 - x0)
  }

  c(meandiff = mean(ediff), maxdiff = max(ediff))
}