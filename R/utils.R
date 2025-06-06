#Function to turn a vector into a string with "," and "and" or "or" for clean messages. 'and.or'
#controls whether words are separated by "and" or "or"; 'is.are' controls whether the list is
#followed by "is" or "are" (to avoid manually figuring out if plural); quotes controls whether
#quotes should be placed around words in string. From WeightIt.
word_list <- function(word.list = NULL, and.or = "and", is.are = FALSE, quotes = FALSE) {
  #When given a vector of strings, creates a string of the form "a and b"
  #or "a, b, and c"
  #If is.are, adds "is" or "are" appropriately

  word.list <- setdiff(word.list, c(NA_character_, ""))

  if (is_null(word.list)) {
    out <- ""
    attr(out, "plural") <- FALSE
    return(out)
  }

  word.list <- add_quotes(word.list, quotes)

  L <- length(word.list)

  if (L == 1L) {
    out <- word.list
    if (is.are) out <- paste(out, "is")
    attr(out, "plural") <- FALSE
    return(out)
  }

  if (is_null(and.or) || isFALSE(and.or)) {
    out <- toString(word.list)
  }
  else {
    and.or <- match_arg(and.or, c("and", "or"))

    if (L == 2L) {
      out <- sprintf("%s %s %s",
                     word.list[1L],
                     and.or,
                     word.list[2L])
    }
    else {
      out <- sprintf("%s, %s %s",
                     toString(word.list[-L]),
                     and.or,
                     word.list[L])
    }
  }

  if (is.are) out <- sprintf("%s are", out)

  attr(out, "plural") <- TRUE

  out
}

#Add quotes to a string
add_quotes <- function(x, quotes = 2L) {
  if (isFALSE(quotes)) {
    return(x)
  }

  if (isTRUE(quotes)) {
    quotes <- '"'
  }

  if (chk::vld_string(quotes)) {
    return(paste0(quotes, x, str_rev(quotes)))
  }

  if (!chk::vld_count(quotes) || quotes > 2L) {
    stop("`quotes` must be boolean, 1, 2, or a string.")
  }

  if (quotes == 0L) {
    return(x)
  }

  x <- {
    if (quotes == 1L) sprintf("'%s'", x)
    else sprintf('"%s"', x)
  }

  x
}

str_rev <- function(x) {
  vapply(lapply(strsplit(x, NULL), rev), paste, character(1L), collapse = "")
}

#More informative and cleaner version of base::match.arg(). Uses chk.
match_arg <- function(arg, choices, several.ok = FALSE) {
  #Replaces match.arg() but gives cleaner error message and processing
  #of arg.
  if (missing(arg)) {
    stop("No argument was supplied to match_arg.")
  }

  arg.name <- deparse1(substitute(arg), width.cutoff = 500L)

  if (missing(choices)) {
    sysP <- sys.parent()
    formal.args <- formals(sys.function(sysP))
    choices <- eval(formal.args[[as.character(substitute(arg))]],
                    envir = sys.frame(sysP))
  }

  if (is_null(arg)) {
    return(choices[1L])
  }

  if (several.ok) {
    chk::chk_character(arg, x_name = add_quotes(arg.name, "`"))
  }
  else {
    chk::chk_string(arg, x_name = add_quotes(arg.name, "`"))
    if (identical(arg, choices)) {
      return(arg[1L])
    }
  }

  i <- pmatch(arg, choices, nomatch = 0L, duplicates.ok = TRUE)
  if (all(i == 0L))
    .err(sprintf("the argument to `%s` should be %s%s",
                 arg.name,
                 ngettext(length(choices), "", if (several.ok) "at least one of " else "one of "),
                 word_list(choices, and.or = "or", quotes = 2L)))
  i <- i[i > 0L]

  choices[i]
}

# Version of interaction(., drop = TRUE) that doesn't succumb to vector limit reached by
# avoiding Cartesian expansion. Falls back to interaction() for small problems.
interaction2 <- function(..., sep = ".", lex.order = TRUE) {

  narg <- ...length()

  if (narg == 0L) {
    stop("No factors specified")
  }

  if (narg == 1L && is.list(..1)) {
    args <- ..1
    narg <- length(args)
  }
  else {
    args <- list(...)
  }

  for (i in seq_len(narg)) {
    args[[i]] <- as.factor(args[[i]])
  }

  if (do.call("prod", lapply(args, nlevels)) <= 1e6) {
    return(interaction(args, drop = TRUE, sep = sep,
                       lex.order = if (is.null(lex.order)) TRUE else lex.order))
  }

  out <- do.call(function(...) paste(..., sep = sep), args)

  args_char <- lapply(args, function(x) {
    x <- unclass(x)
    formatC(x, format = "d", flag = "0", width = ceiling(log10(max(x))))
  })

  lev <- {
    if (is.null(lex.order)) unique(out)
    else if (lex.order) unique(out[order(do.call("paste", c(args_char, list(sep = sep))))])
    else unique(out[order(do.call("paste", c(rev(args_char), list(sep = sep))))])
  }

  factor(out, levels = lev)
}

#Turn a vector into a 0/1 vector. 'zero' and 'one' can be supplied to make it clear which is
#which; otherwise, a guess is used. From WeightIt.
binarize <- function(variable, zero = NULL, one = NULL) {
  var.name <- deparse1(substitute(variable))
  if (is.character(variable) || is.factor(variable)) {
    variable <- factor(variable, nmax = if (is.factor(variable)) nlevels(variable) else NA)
    unique.vals <- levels(variable)
  }
  else {
    unique.vals <- unique(variable)
  }

  if (length(unique.vals) == 1L) {
    return(rep_with(1L, variable))
  }

  if (length(unique.vals) != 2L) {
    .err(sprintf("cannot binarize %s: more than two levels", var.name))
  }

  if (is_not_null(zero)) {
    if (!zero %in% unique.vals) {
      .err(sprintf("the argument to `zero` is not the name of a level of %s", var.name))
    }

    return(setNames(as.integer(variable != zero), names(variable)))
  }

  if (is_not_null(one)) {
    if (!one %in% unique.vals) {
      .err(sprintf("the argument to `one` is not the name of a level of %s", var.name))
    }

    return(setNames(as.integer(variable == one), names(variable)))
  }

  if (is.logical(variable)) {
    return(setNames(as.integer(variable), names(variable)))
  }

  if (is.numeric(variable)) {
    zero <- {
      if (any(unique.vals == 0)) 0
      else min(unique.vals, na.rm = TRUE)
    }

    return(setNames(as.integer(variable != zero), names(variable)))
  }

  variable.numeric <- {
    if (can_str2num(unique.vals)) setNames(str2num(unique.vals), unique.vals)[variable]
    else as.numeric(factor(variable, levels = unique.vals))
  }

  zero <- {
    if (0 %in% variable.numeric) 0
    else min(variable.numeric, na.rm = TRUE)
  }

  setNames(as.integer(variable.numeric != zero), names(variable))
}

is_null <- function(x) length(x) == 0L
is_not_null <- function(x) !is_null(x)

null_or_error <- function(x) {is_null(x) || inherits(x, "try-error")}

#Determine whether a character vector can be coerced to numeric
can_str2num <- function(x) {
  if (is.numeric(x) || is.logical(x)) {
    return(TRUE)
  }

  nas <- is.na(x)
  x_num <- suppressWarnings(as.numeric(as.character(x[!nas])))

  !anyNA(x_num)
}

#Cleanly coerces a character vector to numeric; best to use after can_str2num()
str2num <- function(x) {
  nas <- is.na(x)
  if (!is.numeric(x) && !is.logical(x)) x <- as.character(x)
  x_num <- suppressWarnings(as.numeric(x))
  is.na(x_num)[nas] <- TRUE
  x_num
}

#Capitalize first letter of string
firstup <- function(x) {
  substr(x, 1L, 1L) <- toupper(substr(x, 1L, 1L))
  x
}

#Capitalize first letter of each word
capwords <- function(s, strict = FALSE) {
  cap <- function(s) paste0(toupper(substring(s, 1L, 1L)),
                            {s <- substring(s, 2L)
                            if (strict) tolower(s) else s},
                            collapse = " ")
  sapply(strsplit(s, split = " ", fixed = TRUE), cap,
         USE.NAMES = is_not_null(names(s)))
}

#Reverse a string
str_rev <- function(x) {
  vapply(lapply(strsplit(x, NULL), rev), paste, character(1L), collapse = "")
}

#Clean printing of data frames with numeric and NA elements.
round_df_char <- function(df, digits, pad = "0", na_vals = "") {
  if (NROW(df) == 0L || NCOL(df) == 0L) {
    return(df)
  }

  if (!is.data.frame(df)) {
    df <- as.data.frame.matrix(df, stringsAsFactors = FALSE)
  }

  rn <- rownames(df)
  cn <- colnames(df)

  infs <- o.negs <- array(FALSE, dim = dim(df))
  nas <- is.na(df)
  nums <- vapply(df, is.numeric, logical(1))

  for (i in which(nums)) {
    infs[, i] <- !nas[, i] & !is.finite(df[[i]])
  }

  for (i in which(!nums)) {
    if (can_str2num(df[[i]])) {
      df[[i]] <- str2num(df[[i]])
      nums[i] <- TRUE
    }
  }

  o.negs[, nums] <- !nas[, nums] & df[nums] < 0 & round(df[nums], digits) == 0
  df[nums] <- round(df[nums], digits = digits)

  pad0 <- identical(as.character(pad), "0")

  for (i in which(nums)) {
    df[[i]] <- format(df[[i]], scientific = FALSE, justify = "none", trim = TRUE,
                      drop0trailing = !pad0)

    if (!pad0 && any(grepl(".", df[[i]], fixed = TRUE))) {
      s <- strsplit(df[[i]], ".", fixed = TRUE)
      lengths <- lengths(s)
      digits.r.of.. <- rep.int(0, NROW(df))
      digits.r.of..[lengths > 1] <- nchar(vapply(s[lengths > 1], `[[`, character(1L), 2))

      dots <- rep.int("", length(s))
      dots[lengths <= 1] <- if (as.character(pad) != "") "." else pad

      pads <- vapply(max(digits.r.of..) - digits.r.of..,
                     function(n) paste(rep.int(pad, n), collapse = ""),
                     character(1L))

      df[[i]] <- paste0(df[[i]], dots, pads)
    }
  }

  df[o.negs] <- paste0("-", df[o.negs])

  # Insert NA placeholders
  df[nas] <- na_vals
  df[infs] <- "N/A"

  if (length(rn) > 0) rownames(df) <- rn
  if (length(cn) > 0) names(df) <- cn

  df
}

#Generalized inverse; port of MASS::ginv()
generalized_inverse <- function(sigma, tol = 1e-8) {
  sigmasvd <- svd(sigma)

  pos <- sigmasvd$d > max(tol * sigmasvd$d[1L], 0)

  sigmasvd$v[, pos, drop = FALSE] %*% (sigmasvd$d[pos]^-1 * t(sigmasvd$u[, pos, drop = FALSE]))
}

#(Weighted) variance that uses special formula for binary variables
wvar <- function(x, bin.var = NULL, w = NULL) {
  if (is_null(w)) w <- rep.int(1, length(x))
  if (is_null(bin.var)) bin.var <- all(x == 0 | x == 1)

  w <- w / sum(w) #weights normalized to sum to 1
  mx <- sum(w * x) #weighted mean

  if (bin.var) {
    return(mx * (1 - mx))
  }

  #Reliability weights variance; same as cov.wt()
  sum(w * (x - mx)^2) / (1 - sum(w^2))
}

#Weighted mean faster than weighted.mean()
wm <- function(x, w = NULL, na.rm = TRUE) {
  if (is_null(w)) {
    if (anyNA(x)) {
      if (!na.rm) return(NA_real_)
      nas <- which(is.na(x))
      x <- x[-nas]
    }
    return(sum(x) / length(x))
  }

  if (anyNA(x) || anyNA(w)) {
    if (!na.rm) return(NA_real_)
    nas <- which(is.na(x) | is.na(w))
    x <- x[-nas]
    w <- w[-nas]
  }

  sum(x * w) / sum(w)
}

#Faster diff()
diff1 <- function(x) {
  x[-1L] - x[-length(x)]
}

#cumsum() for probabilities to ensure they are between 0 and 1
.cumsum_prob <- function(x) {
  s <- cumsum(x)
  s / s[length(s)]
}

#Make vector sum to 1, optionally by group
.make_sum_to_1 <- function(x, by = NULL) {
  if (is_null(by)) {
    return(x / sum(x))
  }

  for (i in unique(by)) {
    in_i <- which(by == i)
    x[in_i] <- x[in_i] / sum(x[in_i])
  }

  x
}

#Make vector sum to n (average of 1), optionally by group
.make_sum_to_n <- function(x, by = NULL) {
  if (is_null(by)) {
    return(length(x) * x / sum(x))
  }

  for (i in unique(by)) {
    in_i <- which(by == i)
    x[in_i] <- length(in_i) * x[in_i] / sum(x[in_i])
  }

  x
}

#A faster na.omit for vectors
na.rem <- function(x) {
  x[!is.na(x)]
}

#Extract variables from ..., similar to ...elt() or get0(), by name without evaluating list(...)
...get <- function(x, ifnotfound = NULL) {
  expr <- quote({
    .m1 <- match(.x, ...names())
    if (anyNA(.m1)) {
      .ifnotfound
    }
    else {
      .m2 <- ...elt(.m1[1L])
      if (is_not_null(.m2)) .m2
      else .ifnotfound
    }
  })

  eval(expr,
       pairlist(.x = x[1L], .ifnotfound = ifnotfound),
       parent.frame(1L))
}

#Extract multiple variables from ..., similar to mget(), by name without evaluating list(...)
...mget <- function(x) {
  found <- match(x, eval(quote(...names()), parent.frame(1L)))

  not_found <- is.na(found)

  if (all(not_found)) {
    return(list())
  }

  setNames(lapply(found[!not_found], function(z) {
    eval(quote(...elt(.z)),
         pairlist(.z = z),
         parent.frame(3L))
  }), x[!not_found])
}

#Helper function to fill named vectors with x and given names of y
rep_with <- function(x, y) {
  setNames(rep.int(x, length(y)), names(y))
}

#cat() if verbose = TRUE (default sep = "", line wrapping)
.cat_verbose <- function(..., verbose = TRUE, sep = "") {
  if (!verbose) {
    return(invisible(NULL))
  }

  m <- do.call(function(...) paste(..., sep = sep), list(...))

  if (endsWith(m, "\n")) {
    m <- paste0(paste(strwrap(m), collapse = "\n"), "\n")
  }
  else {
    m <- paste(strwrap(m), collapse = "\n")
  }

  cat(paste(m, collapse = "\n"))
}

#Functions for error handling; based on chk and rlang
pkg_caller_call <- function() {
  pn <- utils::packageName()
  package.funs <- c(getNamespaceExports(pn),
                    .getNamespaceInfo(asNamespace(pn), "S3methods")[, 3L])

  for (i in seq_len(sys.nframe())) {
    e <- sys.call(i)

    n <- rlang::call_name(e)

    if (is_null(n)) {
      next
    }

    if (n %in% package.funs) {
      return(e)
    }
  }

  NULL
}

.err <- function(..., n = NULL, tidy = TRUE) {
  m <- chk::message_chk(..., n = n, tidy = tidy)
  rlang::abort(paste(strwrap(m), collapse = "\n"),
               call = pkg_caller_call())
}
.wrn <- function(..., n = NULL, tidy = TRUE, immediate = TRUE) {
  m <- chk::message_chk(..., n = n, tidy = tidy)

  if (immediate && isTRUE(all.equal(0, getOption("warn")))) {
    rlang::with_options({
      rlang::warn(paste(strwrap(m), collapse = "\n"))
    }, warn = 1)
  }
  else {
    rlang::warn(paste(strwrap(m), collapse = "\n"))
  }
}
.msg <- function(..., n = NULL, tidy = TRUE) {
  m <- chk::message_chk(..., n = n, tidy = tidy)
  rlang::inform(paste(strwrap(m), collapse = "\n"), tidy = FALSE)
}