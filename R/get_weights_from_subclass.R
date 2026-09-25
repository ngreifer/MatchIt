get_weights_from_subclass <- function(subclass, treat, estimand = "ATT", s.weights = NULL) {

  if (is_null(s.weights)) {
    s.weights <- rep_with(1, treat)
  }

  NAsub <- is.na(subclass)

  i1 <- treat == 1 & !NAsub
  i0 <- treat == 0 & !NAsub

  if (sum(s.weights[i1]) == 0) {
    if (sum(s.weights[i0]) == 0) {
      arg::err("no units were matched")
    }

    arg::err("no treated units were matched")
  }
  else if (sum(s.weights[i0]) == 0) {
    arg::err("no control units were matched")
  }

  w <- rep_with(0, treat)

  if (!is.factor(subclass)) {
    subclass <- factor(subclass, nmax = min(sum(i1), sum(i0)))
  }

  #Total sampling weight of the control and treated units in each subclass, in the
  #order of `levels(subclass)`. `split()` keeps the units of each subclass in their
  #original order, so each total is exactly the `sum()` over that subclass alone, but
  #the sample is passed over once rather than once per subclass. Names are dropped
  #because they only slow `split()` down.
  s.weights <- unname(s.weights)

  mass0 <- vapply(split(s.weights[i0], subclass[i0]), sum, numeric(1L))
  mass1 <- vapply(split(s.weights[i1], subclass[i1]), sum, numeric(1L))

  subclass <- unclass(subclass)

  if (estimand == "ATT") {
    w[i1] <- 1
    w[i0] <- (mass1 / mass0)[subclass[i0]]
  }
  else if (estimand == "ATC") {
    w[i1] <- (mass0 / mass1)[subclass[i1]]
    w[i0] <- 1
  }
  else if (estimand == "ATE") {
    w[i1] <- 1 + (mass0 / mass1)[subclass[i1]]
    w[i0] <- 1 + (mass1 / mass0)[subclass[i0]]
  }

  w
}