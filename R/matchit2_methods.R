# this function takes inputs from matchit() and returns the
# strata for each observation in the subclass entry, the
# weight for each observation in the weight entry, and the
# match.matrix object
#
# MATCHIT method = NULL--------------------------------------
matchit2null <- function(discarded, ...) {

  res <- list(weights = as.numeric(!discarded))
  class(res) <- "matchit"

  return(res)
}
# MATCHIT method = cem--------------------------------------
matchit2cem <- function(treat, covs, estimand = "ATT", verbose = FALSE, ...) {
  if (length(covs) == 0) stop("Covariates must be specified in the input formula to use coarsened exact matching.", call. = FALSE)

  if (verbose) cat("Coarsened exact matching...\n")

  A <- list(...)

  estimand <- toupper(estimand)
  estimand <- match_arg(estimand, c("ATT", "ATC", "ATE"))

  #Uses in-house cem, no need for cem package. See cem_matchit.R for code.
  strat <- do.call("cem_matchit", c(list(treat = treat, X = covs, estimand = estimand),
                                    A[names(A) %in% names(formals(cem_matchit))]),
                   quote = TRUE)

  levels(strat) <- seq_len(nlevels(strat))
  names(strat) <- names(treat)

  mm <- NULL
  if (isTRUE(A[["k2k"]])) {
    mm <- nummm2charmm(subclass2mmC(strat, treat, focal = switch(estimand, "ATC" = 0, 1)),
                       treat)
  }

  if (verbose) cat("Calculating matching weights... ")

  res <- list(match.matrix = mm,
              subclass = strat,
              weights = weights.subclass(strat, treat, estimand))

  if (verbose) cat("Done.\n")

  class(res) <- "matchit"
  return(res)
}
# MATCHIT method = exact------------------------------------
matchit2exact <- function(treat, covs, data, estimand = "ATT", verbose = FALSE, ...){

  if(verbose)
    cat("Exact matching... \n")

  if (length(covs) == 0) stop("Covariates must be specified in the input formula to use exact matching.", call. = FALSE)

  estimand <- toupper(estimand)
  estimand <- match_arg(estimand, c("ATT", "ATC", "ATE"))

  xx <- exactify(covs, names(treat))
  cc <- intersect(xx[treat==1], xx[treat==0])

  if (length(cc) == 0) {
    stop("No exact matches were found.", call. = FALSE)
  }

  psclass <- setNames(factor(match(xx, cc), nmax = length(cc)), names(treat))

  if (verbose) cat("Calculating matching weights... ")

  res <- list(subclass = psclass,
              weights = weights.subclass(psclass, treat, estimand))

  if (verbose) cat("Done.\n")

  class(res) <- "matchit"
  return(res)
}

# MATCHIT method = full-------------------------------------
matchit2full <- function(treat, formula, data, distance, discarded,
                         caliper = NULL, mahvars = NULL, exact = NULL,
                         estimand = "ATT", verbose = FALSE,
                         is.full.mahalanobis, min.controls = 0,
                         max.controls = Inf, ...) {

  check.package("optmatch")

  if (verbose) cat("Full matching... \n")

  A <- list(...)

  fm.args <- c("omit.fraction", "mean.controls", "tol")
  A[!names(A) %in% fm.args] <- NULL

  #Set max problem size to Inf and return to original value after match
  omps <- getOption("optmatch_max_problem_size")
  on.exit(options(optmatch_max_problem_size = omps))
  options(optmatch_max_problem_size = Inf)

  estimand <- toupper(estimand)
  estimand <- match_arg(estimand, c("ATT", "ATC", "ATE"))
  if (estimand == "ATC") {
    tc <- c("control", "treated")
    focal <- 0
  }
  else {
    tc <- c("treated", "control")
    focal <- 1
  }

  treat_ <- setNames(as.integer(treat[!discarded] == focal), names(treat)[!discarded])

  if (!is.null(data)) data <- data[!discarded,]

  if (is.full.mahalanobis) {
    if (length(attr(terms(formula, data = data), "term.labels")) == 0) {
      stop("Covariates must be specified in the input formula when distance = \"mahalanobis\".", call. = FALSE)
    }
    mahvars <- formula
  }

  #Create distance matrix; note that Mahalanobis distance computed using entire
  #sample, like method2nearest, as opposed to within exact strata, like optmatch.
  if (!is.null(mahvars)) {
    mahvars <- update(mahvars, treat_ ~ .)
    environment(mahvars) <- sys.frame(sys.nframe())

    mo <- optmatch::match_on(mahvars, data = data[names(data) != "treat_"],
                             method = "mahalanobis")
  }
  else {
    mo <- optmatch::match_on(distance[!discarded], z = treat_)
  }

  if (!is.null(exact)) {
    ex <- exactify(model.frame(exact, data = data))

    cc <- intersect(ex[treat_==1], ex[treat_==0])
    if (length(cc) == 0) stop("No matches were found.", call. = FALSE)

    mo <- mo + optmatch::exactMatch(ex, treat_)
  }

  if (!is.null(caliper)) {
    if (min.controls != 0) {
      stop("Calipers cannot be used with method = \"full\" when 'min.controls' is specified.", call. = FALSE)
    }

    if (any(names(caliper) != "")) {
      cov.cals <- setdiff(names(caliper), "")
      calcovs <- get.covs.matrix(reformulate(cov.cals, intercept = FALSE), data = data)
    }
    for (i in seq_along(caliper)) {
      if (names(caliper)[i] == "") {
        mo_cal <- optmatch::match_on(distance[!discarded], z = treat_)
      }
      else {
        mo_cal <- optmatch::match_on(setNames(calcovs[,names(caliper)[i]], names(treat_)), z = treat_)
      }

      mo <- mo + optmatch::caliper(mo_cal, caliper[i])
    }
    rm(mo_cal)
  }

  pair <- setNames(rep(NA_character_, length(treat)), names(treat))

  withCallingHandlers({
    p <- do.call(optmatch::fullmatch,
                 c(list(mo,
                        min.controls = min.controls,
                        max.controls = max.controls,
                        data = treat), #just to get rownames; not actually used in matching
                   A))
  },
  warning = function(w) {
    warning(paste0("(from optmatch) ", conditionMessage(w)), call. = FALSE, immediate. = TRUE)
    invokeRestart("muffleWarning")
  },
  error = function(e) {
    stop(paste0("(from optmatch) ", conditionMessage(e)), call. = FALSE)
  })

  pair[names(p)[!is.na(p)]] <- as.character(p[!is.na(p)])

  if (all(is.na(pair))) stop("No matches were found.", call. = FALSE)

  psclass <- factor(pair, nmax = min(sum(treat_ == 1), sum(treat_ == 0)))
  levels(psclass) <- seq_len(nlevels(psclass))
  names(psclass) <- names(treat)

  #No match.matrix because treated units don't index matched strata (i.e., more than one
  #treated unit can be in the same stratum). Stratum information is contained in subclass.

  if (verbose) cat("Calculating matching weights... ")

  res <- list(subclass = psclass,
              weights = weights.subclass(psclass, treat, estimand))

  if (verbose) cat("Done.\n")

  class(res) <- c("matchit")
  return(res)
}

# MATCHIT method = optimal----------------------------------
matchit2optimal <- function(treat, formula, data, distance, discarded,
                            ratio = 1, caliper = NULL, mahvars = NULL, exact = NULL,
                            estimand = "ATT", verbose = FALSE,
                            is.full.mahalanobis, min.controls = NULL,
                            max.controls = NULL, ...) {

  check.package("optmatch")

  if (verbose) cat("Optimal matching... \n")

  A <- list(...)
  pm.args <- c("tol")
  A[!names(A) %in% pm.args] <- NULL

  #Set max problem size to Inf and return to original value after match
  omps <- getOption("optmatch_max_problem_size")
  on.exit(options(optmatch_max_problem_size = omps))
  options(optmatch_max_problem_size = Inf)

  r <- process.ratio(ratio, min.controls, max.controls)
  ratio <- r["ratio"]
  if (!is.null(max.controls)) {
    min.controls <- r["min.controls"]
    max.controls <- r["max.controls"]
  }

  estimand <- toupper(estimand)
  estimand <- match_arg(estimand, c("ATT", "ATC"))
  if (estimand == "ATC") {
    tc <- c("control", "treated")
    focal <- 0
  }
  else {
    tc <- c("treated", "control")
    focal <- 1
  }

  treat_ <- setNames(as.integer(treat[!discarded] == focal), names(treat)[!discarded])

  if (!is.null(data)) data <- data[!discarded,]

  if (is.full.mahalanobis) {
    if (length(attr(terms(formula, data = data), "term.labels")) == 0) {
      stop("Covariates must be specified in the input formula when distance = \"mahalanobis\".", call. = FALSE)
    }
    mahvars <- formula
  }

  if (!is.null(caliper)) {
    warning("Calipers are currently not compatible with method = \"optimal\" and will be ignored.", call. = FALSE, immediate. = TRUE)
    caliper <- NULL
  }

  if (is.null(max.controls)) {
    min.controls <- max.controls <- ratio
  }

  #Exact matching strata
  if (!is.null(exact)) {
    ex <- factor(exactify(model.frame(exact, data = data)))

    cc <- intersect(ex[treat_==1], ex[treat_==0])
    if (length(cc) == 0) stop("No matches were found.", call. = FALSE)

    e_ratios <- sapply(levels(ex), function(e) sum(treat_[ex == e] == 0)/sum(treat_[ex == e] == 1))

    if (any(e_ratios < 1)) {
      warning(paste0("Fewer ", tc[2], " units than ", tc[1], " units in some 'exact' strata; not all ", tc[1], " units will get a match."), immediate. = TRUE, call. = FALSE)
    }
    if (ratio > 1 && any(e_ratios < ratio)) {
      if (ratio == max.controls)
        warning(paste0("Not all ", tc[1], " units will get ", ratio, " matches."), immediate. = TRUE, call. = FALSE)
      else
        warning(paste0("Not enough ", tc[2], " units for an average of ", ratio, " matches per ", tc[1], " unit in all 'exact' strata."), immediate. = TRUE, call. = FALSE)
    }
  }
  else {
    ex <- factor(rep("_", length(treat_)), levels = "_")
    e_ratios <- setNames(sum(treat_ == 0)/sum(treat_ == 1), levels(ex))

    if (e_ratios < 1) {
      warning(paste0("Fewer ", tc[2], " units than ", tc[1], " units; not all ", tc[1], " units will get a match."), immediate. = TRUE, call. = FALSE)
    }
    else if (e_ratios < ratio) {
      if (ratio == max.controls)
        warning(paste0("Not all ", tc[1], " units will get ", ratio, " matches."), immediate. = TRUE, call. = FALSE)
      else
        warning(paste0("Not enough ", tc[2], " units for an average of ", ratio, " matches per ", tc[1], " unit."), immediate. = TRUE, call. = FALSE)
    }
  }

  #Create distance matrix; note that Mahalanobis distance computed using entire
  #sample, like method2nearest, as opposed to within exact strata, like optmatch.
  if (!is.null(mahvars)) {
    mahvars <- update(mahvars, treat_ ~ .)
    environment(mahvars) <- sys.frame(sys.nframe())

    mo <- optmatch::match_on(mahvars, data = data[names(data) != "treat_"],
                             method = "mahalanobis")
  }
  else {
    mo <- optmatch::match_on(distance[!discarded], z = treat_)
  }

  #Initialize pair membership; must include names
  pair <- setNames(rep(NA_character_, length(treat)), names(treat))

  for (e in levels(ex)) {
    if (nlevels(ex) > 1) mo_ <- mo[ex[treat_==1] == e, ex[treat_==0] == e, drop = FALSE]

    #Process ratio, etc., when available ratio in exact matching categories
    #(e_ratio) differs from requested ratio
    if (e_ratios[e] < 1) {
      #Switch treatment and control labels; unmatched treated units are dropped
      ratio_ <- min.controls_ <- max.controls_ <- 1
      mo_ <- t(mo_)
    }
    else if (e_ratios[e] < ratio) {
      #Lower ratio and min.controls.
      ratio_ <- e_ratios[e]
      min.controls_ <- min(min.controls, floor(e_ratios[e]))
      max.controls_ <- max.controls
    }
    else {
      ratio_ <- ratio
      min.controls_ <- min.controls
      max.controls_ <- max.controls
    }

    withCallingHandlers({
        p <- do.call(optmatch::fullmatch,
                        c(list(if (nlevels(ex) > 1) mo_ else mo,
                               mean.controls = ratio_,
                               min.controls = min.controls_,
                               max.controls = max.controls_,
                               data = treat), #just to get rownames; not actually used in matching
                          A))
    },
    warning = function(w) {
      warning(paste0("(from optmatch) ", conditionMessage(w)), call. = FALSE, immediate. = TRUE)
      invokeRestart("muffleWarning")
    },
    error = function(e) {
      stop(paste0("(from optmatch) ", conditionMessage(e)), call. = FALSE)
    })

    pair[names(p)[!is.na(p)]] <- paste(as.character(p[!is.na(p)]), e, sep = "|")
  }

  if (all(is.na(pair))) stop("No matches were found.", call. = FALSE)

  psclass <- factor(pair); levels(psclass) <- seq_len(nlevels(psclass))
  names(psclass) <- names(treat)

  mm <- nummm2charmm(subclass2mmC(psclass, treat, focal), treat)

  if (verbose) cat("Calculating matching weights... ")

  ## calculate weights and return the results
  res <- list(match.matrix = mm,
              subclass = psclass,
              weights = weights.subclass(psclass, treat, estimand))

  if (verbose) cat("Done.\n")

  class(res) <- "matchit"
  return(res)
}
# MATCHIT method = genetic----------------------------------
matchit2genetic <- function(treat, data, distance, discarded,
                            ratio = 1, s.weights = NULL, replace = FALSE, m.order = NULL,
                            caliper = NULL, mahvars = NULL, exact = NULL,
                            formula = NULL, estimand = "ATT", verbose = FALSE,
                            is.full.mahalanobis, use.genetic = TRUE, ...) {

  check.package(c("Matching", "rgenoud"))

  if (verbose) cat("Genetic matching... \n")

  A <- list(...)

  estimand <- toupper(estimand)
  estimand <- match_arg(estimand, c("ATT", "ATC"))
  if (estimand == "ATC") {
    tc <- c("control", "treated")
    focal <- 0
  }
  else {
    tc <- c("treated", "control")
    focal <- 1
  }

  if (!replace) {
    if (sum(!discarded & treat != focal) < sum(!discarded & treat == focal)) {
      warning("Fewer ", tc[2], " units than ", tc[1], " units; not all ", tc[1], " units will get a match.", immediate. = TRUE, call. = FALSE)    }
    else if (sum(!discarded & treat != focal) < sum(!discarded & treat == focal)*ratio) {
      stop(paste0("Not enough ", tc[2], " units for ", ratio, " matches for each ", tc[1], " unit."), call. = FALSE)
    }
  }

  treat <- setNames(as.integer(treat == focal), names(treat))

  ratio <- process.ratio(ratio)

  n.obs <- length(treat)
  n1 <- sum(treat == 1)

  if (is.null(names(treat))) names(treat) <- seq_len(n.obs)

  if (!is.null(distance)) {
    if (is.null(m.order)) m.order <- if (estimand == "ATC") "smallest" else "largest"
    else m.order <- match_arg(m.order, c("largest", "smallest", "random", "data"))
    ord <- switch(m.order,
                  "largest" = order(distance, decreasing = TRUE),
                  "smallest" = order(distance),
                  "random" = sample(seq_len(n.obs), n.obs, replace = FALSE),
                  "data" = seq_len(n.obs))
  }
  else {
    m.order <- match_arg(m.order, c("data", "random"))
    ord <- switch(m.order,
                  "random" = sample(seq_len(n.obs), n.obs, replace = FALSE),
                  "data" = seq_len(n.obs))
  }
  ord <- ord[!ord %in% which(discarded)]

  #Create X (matching variables) and covs_to_balance
  if (!is.null(mahvars)) {
    covs_to_balance <- get.covs.matrix(formula, data = data)
    X <- get.covs.matrix(mahvars, data = data)
  }
  else if (is.full.mahalanobis) {
    covs_to_balance <- get.covs.matrix(formula, data = data)
    X <- covs_to_balance
  }
  else {
    covs_to_balance <- get.covs.matrix(formula, data = data)
    X <- cbind(covs_to_balance, distance)
  }

  if (ncol(covs_to_balance) == 0) {
    stop("Covariates must be specified in the input formula to use genetic matching.", call. = FALSE)
  }

  #Process exact; exact.log will be supplied to GenMatch() and Match()
  if (!is.null(exact)) {
    #Add covariates in exact not in X to X
    ex <- as.integer(factor(exactify(model.frame(exact, data = data), names(treat))))

    cc <- intersect(ex[treat==1], ex[treat==0])
    if (length(cc) == 0) stop("No matches were found.", call. = FALSE)

    X <- cbind(X, ex)

    exact.log <- c(rep(FALSE, ncol(X) - 1), TRUE)
  }
  else exact.log <- ex <- NULL

  #Process caliper; cal will be supplied to GenMatch() and Match()
  if (!is.null(caliper)) {
    #Add covariates in caliper other than distance (cov.cals) not in X to X
    cov.cals <- setdiff(names(caliper), "")
    if (length(cov.cals) > 0 && any(!cov.cals %in% colnames(X))) {
      calcovs <- get.covs.matrix(reformulate(cov.cals[!cov.cals %in% colnames(X)]), data = data)
      X <- cbind(X, calcovs)

      #Expand exact.log for newly added covariates
      if (!is.null(exact.log)) exact.log <- c(exact.log, rep(FALSE, ncol(calcovs)))
    }
    else cov.cals <- NULL

    #Matching::Match multiplies calipers by pop SD, so we need to divide by pop SD to unstandardize
    pop.sd <- function(x) sqrt(sum((x-mean(x))^2)/length(x))
    caliper <- caliper / vapply(names(caliper), function(x) {
      if (x == "") pop.sd(distance[!discarded])
      else pop.sd(X[!discarded, x])
    }, numeric(1L))

    #cal needs one value per variable in X
    cal <- setNames(rep(Inf, ncol(X)), colnames(X))

    #First put covariate calipers into cal
    if (length(cov.cals) > 0) {
      cal[intersect(cov.cals, names(cal))] <- caliper[intersect(cov.cals, names(cal))]
    }

    #Then put distance caliper into cal
    if ("" %in% names(caliper)) {
      dist.cal <- caliper[names(caliper) == ""]
      if (!is.null(mahvars)) {
        #If mahvars specified, distance is not yet in X, so add it to X
        X <- cbind(X, distance)
        cal <- c(cal, dist.cal)

        #Expand exact.log for newly added distance
        if (!is.null(exact.log)) exact.log <- c(exact.log, FALSE)
      }
      else {
        #Otherwise, distance is in X at the specified index
        cal[ncol(covs_to_balance) + 1] <- dist.cal
      }
    }
    else dist.cal <- NULL

  }
  else {
    cal <- dist.cal <- cov.cals <- NULL
  }

  #Reorder data according to m.order since Match matches in order of data;
  #ord already excludes discarded units

  treat_ <- treat[ord]
  covs_to_balance <- covs_to_balance[ord,,drop = FALSE]
  X_ <- X[ord,,drop = FALSE]
  if (!is.null(s.weights)) s.weights <- s.weights[ord]

  if (use.genetic) {
    withCallingHandlers({
      g.out <- do.call(Matching::GenMatch,
                       c(list(Tr = treat_, X = X_, BalanceMatrix = covs_to_balance,
                              M = ratio, exact = exact.log, caliper = cal,
                              replace = replace, estimand = "ATT", ties = FALSE,
                              CommonSupport = FALSE, verbose = verbose,
                              weights = s.weights, print.level = 2*verbose),
                         A[names(A) %in% names(formals(Matching::GenMatch))]))
    },
    warning = function(w) {
      if (!startsWith(conditionMessage(w), "replace==FALSE, but there are more (weighted) treated obs than control obs.")) {
        warning(paste0("(from Matching) ", conditionMessage(w)), call. = FALSE, immediate. = TRUE)
      }
      invokeRestart("muffleWarning")
    },
    error = function(e) {
      stop(paste0("(from Matching) ", conditionMessage(e)), call. = FALSE)
    })
  }
  else {
    #For debugging
    g.out <- NULL
  }

  lab <- names(treat)
  lab1 <- names(treat[treat == 1])

  lab_ <- names(treat_)

  ind_ <- seq_along(treat)[ord]

  # if (!isFALSE(A$use.Match)) {
  withCallingHandlers({
    m.out <- Matching::Match(Tr = treat_, X = X_,
                             M = ratio, exact = exact.log, caliper = cal,
                             replace = replace, estimand = "ATT", ties = FALSE,
                             weights = s.weights, CommonSupport = FALSE, Weight = 3,
                             Weight.matrix = if (use.genetic) g.out
                             else if (is.null(s.weights)) generalized_inverse(cor(X_))
                             else generalized_inverse(cov.wt(X_, s.weights, cor = TRUE)$cor),
                             version = "fast")
  },
  warning = function(w) {
    if (!startsWith(conditionMessage(w), "replace==FALSE, but there are more (weighted) treated obs than control obs.")) {
      warning(paste0("(from Matching) ", conditionMessage(w)), call. = FALSE, immediate. = TRUE)
    }
    invokeRestart("muffleWarning")
  },
  error = function(e) {
    stop(paste0("(from Matching) ", conditionMessage(e)), call. = FALSE)
  })

  #Note: must use character match.matrix because of re-ordering treat into treat_
  mm <- matrix(NA_integer_, nrow = n1, ncol = max(table(m.out$index.treated)),
               dimnames = list(lab1, NULL))

  unique.matched.focal <- unique(m.out$index.treated, nmax = n1)

  ind1__ <- match(lab_, lab1)
  for (i in unique.matched.focal) {
    matched.units <- ind_[m.out$index.control[m.out$index.treated == i]]
    mm[ind1__[i], seq_along(matched.units)] <- matched.units
  }

  # }
  # else {
  #   ord1 <- ord[ord %in% which(treat == 1)]
  #   if (!is.null(cov.cals)) calcovs <- get.covs.matrix(reformulate(cov.cals), data = data)
  #   else calcovs <- NULL
  #
  #   if (is.null(g.out)) MWM <- generalized_inverse(cov(X))
  #   else MWM <- g.out$Weight.matrix %*% diag(1/apply(X, 2, var))
  #
  #   if (isFALSE(A$fast)) {
  #     mm <- nn_match(treat, ord1, ratio, replace, discarded, distance, ex, dist.cal,
  #                    cov.cals, calcovs, X, MWM)
  #   }
  #   else {
  #     mm <- nn_matchC(treat, ord1, ratio, replace, discarded, distance, ex, dist.cal,
  #                     cov.cals, calcovs, X, MWM)
  #   }
  #
  #   mm[] <- names(treat)[mm]
  #   dimnames(mm) <- list(lab1, seq_len(ratio))
  # }

  if (verbose) cat("Calculating matching weights... ")

  if (replace) {
    psclass <- NULL
  }
  else {
    psclass <- mm2subclass(mm, treat)
  }

  if (verbose) cat("Done.\n")

  res <- list(match.matrix = nummm2charmm(mm, treat),
              subclass = psclass,
              weights = weights.matrix(mm, treat),
              obj = g.out)

  class(res) <- "matchit"
  return(res)
}

# MATCHIT method = nearest----------------------------------
matchit2nearest <-  function(treat, data, distance, discarded,
                             ratio = 1, s.weights = NULL, replace = FALSE, m.order = NULL,
                             caliper = NULL, mahvars = NULL, exact = NULL,
                             formula = NULL, estimand = "ATT", verbose = FALSE,
                             is.full.mahalanobis, fast = TRUE,
                             min.controls = NULL, max.controls = NULL, ...){

  if (verbose) {
    if (fast) check.package("RcppProgress")
    cat("Nearest neighbor matching... \n")
  }

  r <- process.ratio(ratio, min.controls, max.controls)
  ratio <- r["ratio"]
  if (!is.null(max.controls)) {
    min.controls <- r["min.controls"]
    max.controls <- r["max.controls"]
  }

  estimand <- toupper(estimand)
  estimand <- match_arg(estimand, c("ATT", "ATC"))
  if (estimand == "ATC") {
    tc <- c("control", "treated")
    focal <- 0
  }
  else {
    tc <- c("treated", "control")
    focal <- 1
  }

  treat <- setNames(as.integer(treat == focal), names(treat))

  n.obs <- length(treat)
  n1 <- sum(treat == 1)

  lab <- names(treat)
  lab1 <- lab[treat == 1]

  if (!is.null(distance)) {
    names(distance) <- names(treat)
  }

  if (!is.null(caliper)) {
    if (any(names(caliper) != "")) {
      caliper.covs <- caliper[names(caliper) != ""]
      caliper.covs.mat <- get.covs.matrix(reformulate(names(caliper.covs)), data = data)
    }
    else {
      caliper.covs.mat <- caliper.covs <- NULL
    }

    if (any(names(caliper) == "")) {
      caliper.dist <- caliper[names(caliper) == ""]
    }
    else {
      caliper.dist <- NULL
    }
  }
  else {
    caliper.dist <- caliper.covs <- NULL
    caliper.covs.mat <- NULL
  }

  if (is.full.mahalanobis) {
    mahcovs <- get.covs.matrix(formula, data)
    if (ncol(mahcovs) == 0) stop("Covariates must be specified in the input formula when distance = \"mahalanobis\".", call. = FALSE)
    if (is.null(s.weights))
      mahSigma_inv <- generalized_inverse(cov(mahcovs))
    else
      mahSigma_inv <- generalized_inverse(cov.wt(mahcovs, s.weights)$cov)
  }
  else if (!is.null(mahvars)) {
    mahcovs <- get.covs.matrix(mahvars, data)
    if (is.null(s.weights))
      mahSigma_inv <- generalized_inverse(cov(mahcovs))
    else
      mahSigma_inv <- generalized_inverse(cov.wt(mahcovs, s.weights)$cov)
  }
  else {
    mahcovs <- mahSigma_inv <- NULL
  }

  if (replace) {
    m.order <- "data"
  }

  if (!is.null(exact)) {
    ex <- factor(exactify(model.frame(exact, data = data)))

    cc <- intersect(ex[treat==1], ex[treat==0])
    if (length(cc) == 0) stop("No matches were found.", call. = FALSE)

    e_ratios <- sapply(levels(ex), function(e) sum(treat[ex == e] == 0)/sum(treat[ex == e] == 1))

    if (any(e_ratios < 1)) {
      warning(paste0("Fewer ", tc[2], " units than ", tc[1], " units in some 'exact' strata; not all ", tc[1], " units will get a match."), immediate. = TRUE, call. = FALSE)
    }
    if (ratio > 1 && any(e_ratios < ratio)) {
      if (is.null(max.controls) || ratio == max.controls)
        warning(paste0("Not all ", tc[1], " units will get ", ratio, " matches."), immediate. = TRUE, call. = FALSE)
      else
        warning(paste0("Not enough ", tc[2], " units for an average of ", ratio, " matches per ", tc[1], " unit in all 'exact' strata."), immediate. = TRUE, call. = FALSE)
    }

    ex <- setNames(factor(ex), lab)
  }
  else {
    ex <- NULL
    # ex <- factor(rep("_", length(treat)), levels = "_")
    e_ratios <- setNames(sum(treat == 0)/sum(treat == 1), levels(ex))

    if (e_ratios < 1) {
      warning(paste0("Fewer ", tc[2], " units than ", tc[1], " units; not all ", tc[1], " units will get a match."), immediate. = TRUE, call. = FALSE)
    }
    else if (e_ratios < ratio) {
      if (is.null(max.controls) || ratio == max.controls)
        warning(paste0("Not all ", tc[1], " units will get ", ratio, " matches."), immediate. = TRUE, call. = FALSE)
      else
        warning(paste0("Not enough ", tc[2], " units for an average of ", ratio, " matches per ", tc[1], " unit."), immediate. = TRUE, call. = FALSE)
    }

  }

  #Variable ratio (extremal matching), Ming & Rosenbaum (2000)
  #Each treated unit get its own value of ratio
  if (!is.null(max.controls)) {
    if (is.null(distance)) stop("'distance' cannot be \"mahalanobis\" for variable ratio matching.", call. = FALSE)

    n1 <- sum(treat == 1)
    m <- round(ratio * n1)
    # if (m > sum(treat == 0)) stop("'ratio' must be less than or equal to n0/n1.", call. = FALSE)

    kmax <- floor((m - min.controls*(n1-1)) / (max.controls - min.controls))
    kmin <- n1 - kmax - 1
    kmed <- m - (min.controls*kmin + max.controls*kmax)

    ratio0 <- c(rep(min.controls, kmin), kmed, rep(max.controls, kmax))

    #Make sure no units are assigned 0 matches
    if (any(ratio0 == 0)) {
      ind <- which(ratio0 == 0)
      ratio0[ind] <- 1
      ratio0[ind + 1] <- ratio0[ind + 1] - 1
    }

    ratio <- rep(NA_integer_, n1)

    #Order by distance; treated are supposed to have higher values
    ratio[order(distance[treat == 1],
                decreasing = mean(distance[treat == 1]) > mean(distance[treat != 1]))] <- ratio0
    ratio <- as.integer(ratio)
  }
  else {
    ratio <- as.integer(rep(ratio, n1))
  }
  max_rat <- max(ratio)

  if (!is.null(distance)) {
    if (is.null(m.order)) m.order <- if (estimand == "ATC") "smallest" else "largest"
    else m.order <- match_arg(m.order, c("largest", "smallest", "random", "data"))
  }
  else {
    m.order <- match_arg(m.order, c("data", "random"))
  }

  if (is.null(ex)) {
    ord <- switch(m.order,
                   "largest" = order(distance[treat == 1], decreasing = TRUE),
                   "smallest" = order(distance[treat == 1]),
                   "random" = sample(seq_len(n1), n1, replace = FALSE),
                   "data" = seq_len(n1))

    mm <- nn_matchC(treat, ord, ratio, max_rat, replace, discarded, distance, ex, caliper.dist,
                    caliper.covs, caliper.covs.mat, mahcovs, mahSigma_inv, verbose)
  }
  else {
    mm_list <- lapply(levels(ex), function(e) {
      .e <- which(ex == e)
      treat_ <- treat[.e]
      discarded_ <- discarded[.e]
      distance_ <- distance[.e]
      caliper.covs.mat_ <- caliper.covs.mat[.e,,drop = FALSE]
      mahcovs_ <- mahcovs[.e,,drop = FALSE]
      ratio_ <- ratio[ex[treat==1]==e]

      n1_ <- sum(treat_ == 1)
      ord_ <- switch(m.order,
                     "largest" = order(distance_[treat_ == 1], decreasing = TRUE),
                     "smallest" = order(distance_[treat_ == 1]),
                     "random" = sample(seq_len(n1_), n1_, replace = FALSE),
                     "data" = seq_len(n1_))

      mm_ <- nn_matchC(treat_, ord_, ratio_, max_rat, replace, discarded_, distance_, NULL, caliper.dist,
                       caliper.covs, caliper.covs.mat_, mahcovs_, mahSigma_inv, verbose)

      mm_[] <- seq_along(treat)[.e][mm_]
      mm_
    })

    mm <- do.call("rbind", mm_list)[lab1,, drop = FALSE]
  }

  if (verbose) cat("Calculating matching weights... ")

  if (replace) {
    psclass <- NULL
  }
  else {
    psclass <- mm2subclass(mm, treat)
  }

  if (verbose) cat("Done.\n")

  res <- list(match.matrix = nummm2charmm(mm, treat),
              subclass = psclass,
              weights = weights.matrix(mm, treat))

  class(res) <- "matchit"

  return(res)
}

# MATCHIT method = subclass---------------------------------
matchit2subclass <- function(treat, distance, discarded,
                             replace = FALSE, exact = NULL,
                             estimand = "ATT", verbose = FALSE,
                             ...) {

  if(verbose)
    cat("Subclassifying... \n")

  A <- list(...)
  subclass <- A[["subclass"]]
  sub.by <- A[["sub.by"]]
  min.n <- A[["min.n"]]

  #Checks
  if (is.null(subclass)) subclass <- 6
  else if (!is.numeric(subclass) || !is.null(dim(subclass))) {
    stop("subclass must be a numeric value.", call. = FALSE)
  }
  else if (length(subclass) == 1) {
    if (round(subclass) <= 1) {
      stop("subclass must be greater than 1.",call.=FALSE)
    }
  }
  else if (!all(subclass <= 1 & subclass >= 0)) {
    stop("When specifying subclass as a vector of quantiles, all values must be between 0 and 1.",
         call. = FALSE)
  }

  if (!is.null(sub.by)) {
    sub.by.choices <- c("treat", "control", "all")
    if (!is.character(sub.by) || length(sub.by) != 1 || anyNA(pmatch(sub.by, sub.by.choices))) {
      stop("'sub.by' is deprecated and can't be converted into a proper input. Please supply an argument to 'estimand' instead.", call. = FALSE)
    }
    else {
      sub.by <- sub.by.choices[pmatch(sub.by, sub.by.choices)]
      estimand <- switch(sub.by, "treat" = "ATT", "control" = "ATC", "ATE")
      warning(paste0("'sub.by' is deprecated and has been replaced with 'estimand'. Setting 'estimand' to \"", estimand, "\"."), call. = FALSE, immediate. = TRUE)
    }
  }
  else {
    estimand <- toupper(estimand)
    estimand <- match_arg(estimand, c("ATT", "ATC", "ATE"))
  }

  if (is.null(min.n)) min.n <- 1
  else if (!is.numeric(min.n) || length(min.n) != 1) {
    stop("'min.n' must be a single number.", call. = FALSE)
  }

  n.obs <- length(treat)

  ## Setting Cut Points
  if (length(subclass) == 1) {
    sprobs <- seq(0, 1, length.out = round(subclass) + 1)
  }
  else {
    sprobs <- sort(subclass)
    if (sprobs[1] != 0) sprobs <- c(0, sprobs)
    if (sprobs[length(sprobs)] != 1) sprobs <- c(sprobs, 1)
    subclass <- length(sprobs) - 1
  }

  q <- switch(estimand,
              "ATT" = quantile(distance[treat==1], probs = sprobs, na.rm = TRUE),
              "ATC" = quantile(distance[treat==0], probs = sprobs, na.rm = TRUE),
              "ATE" = quantile(distance, probs = sprobs, na.rm = TRUE))

  ## Calculating Subclasses
  psclass <- setNames(rep(NA_integer_, n.obs), names(treat))
  psclass[!discarded] <- as.integer(findInterval(distance[!discarded], q, all.inside = TRUE))

  ## If any subclasses don't have members of a treatment group, fill them
  ## by "scooting" units from nearby subclasses until each subclass has a unit
  ## from each treatment group
  if (any(table(treat, psclass) < min.n)) {
    psclass[!discarded] <- subclass_scoot(psclass[!discarded], treat[!discarded], distance[!discarded], min.n)
  }

  psclass <- setNames(factor(psclass, nmax = length(q)), names(treat))

  if (verbose) cat("Done\n")

  #warning for discrete data

  if (nlevels(psclass) != subclass){
    warning("Due to discreteness in the distance measure, fewer subclasses were generated than were requested.", call.=FALSE)
  }

  if (verbose) cat("Calculating matching weights... ")

  res <- list(subclass = psclass, q.cut = q,
              qn = qn(treat, psclass, discarded),
              weights = weights.subclass(psclass, treat, estimand))

  if (verbose) cat("Done.\n")

  class(res) <- c("matchit.subclass", "matchit")
  return(res)
}

