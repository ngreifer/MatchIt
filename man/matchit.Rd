% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/matchit.R
\name{matchit}
\alias{matchit}
\title{Matching for Causal Inference}
\usage{
matchit(
  formula,
  data = NULL,
  method = "nearest",
  distance = "glm",
  link = "logit",
  distance.options = list(),
  estimand = "ATT",
  exact = NULL,
  mahvars = NULL,
  antiexact = NULL,
  discard = "none",
  reestimate = FALSE,
  s.weights = NULL,
  replace = FALSE,
  m.order = NULL,
  caliper = NULL,
  std.caliper = TRUE,
  ratio = 1,
  verbose = FALSE,
  include.obj = FALSE,
  normalize = TRUE,
  ...
)
}
\arguments{
\item{formula}{a two-sided \code{\link{formula}} object containing the treatment and
covariates to be used in creating the distance measure used in the matching.
This formula will be supplied to the functions that estimate the distance
measure. The formula should be specified as \code{A ~ X1 + X2 + ...} where
\code{A} represents the treatment variable and \code{X1} and \code{X2} are
covariates.}

\item{data}{a data frame containing the variables named in \code{formula}
and possible other arguments. If not found in \code{data}, the variables
will be sought in the environment.}

\item{method}{the matching method to be used. The allowed methods are
\code{\link[=method_nearest]{"nearest"}} for nearest neighbor matching (on
the propensity score by default), \code{\link[=method_optimal]{"optimal"}}
for optimal pair matching, \code{\link[=method_full]{"full"}} for optimal
full matching, \code{\link[=method_quick]{"quick"}} for generalized (quick)
full matching, \code{\link[=method_genetic]{"genetic"}} for genetic
matching, \code{\link[=method_cem]{"cem"}} for coarsened exact matching,
\code{\link[=method_exact]{"exact"}} for exact matching,
\code{\link[=method_cardinality]{"cardinality"}} for cardinality and
profile matching, and \code{\link[=method_subclass]{"subclass"}} for
subclassification. When set to \code{NULL}, no matching will occur, but
propensity score estimation and common support restrictions will still occur
if requested. See the linked pages for each method for more details on what
these methods do, how the arguments below are used by each on, and what
additional arguments are allowed.}

\item{distance}{the distance measure to be used. Can be either the name of a
method of estimating propensity scores (e.g., \code{"glm"}), the name of a
method of computing a distance matrix from the covariates (e.g.,
\code{"mahalanobis"}), a vector of already-computed distance measures, or a
matrix of pairwise distances. See \code{\link{distance}} for allowable
options. The default is \code{"glm"} for propensity scores estimated with
logistic regression using \code{\link[=glm]{glm()}}. Ignored for some methods; see individual
methods pages for information on whether and how the distance measure is
used.}

\item{link}{when \code{distance} is specified as a string, an additional
argument controlling the link function used in estimating the distance
measure. Allowable options depend on the specific \code{distance} value
specified. See \code{\link{distance}} for allowable options with each
option. The default is \code{"logit"}, which, along with \code{distance = "glm"}, identifies the default measure as logistic regression propensity scores.}

\item{distance.options}{a named list containing additional arguments
supplied to the function that estimates the distance measure as determined
by the argument to \code{distance}. See \code{\link{distance}} for an
example of its use.}

\item{estimand}{a string containing the name of the target estimand desired.
Can be one of \code{"ATT"}, \code{"ATC"}, or \code{"ATE"}. Default is \code{"ATT"}. See Details and the individual methods
pages for information on how this argument is used.}

\item{exact}{for methods that allow it, for which variables exact matching
should take place. Can be specified as a string containing the names of
variables in \code{data} to be used or a one-sided formula with the desired
variables on the right-hand side (e.g., \code{~ X3 + X4}). See the
individual methods pages for information on whether and how this argument is
used.}

\item{mahvars}{for methods that allow it, on which variables Mahalanobis
distance matching should take place when \code{distance} corresponds to
propensity scores. Usually used to perform Mahalanobis distance matching
within propensity score calipers, where the propensity scores are computed
using \code{formula} and \code{distance}. Can be specified as a string
containing the names of variables in \code{data} to be used or a one-sided
formula with the desired variables on the right-hand side (e.g., \code{~ X3 + X4}). See the individual methods pages for information on whether and how this argument is used.}

\item{antiexact}{for methods that allow it, for which variables anti-exact
matching should take place. Anti-exact matching ensures paired individuals
do not have the same value of the anti-exact matching variable(s). Can be
specified as a string containing the names of variables in \code{data} to be
used or a one-sided formula with the desired variables on the right-hand
side (e.g., \code{~ X3 + X4}). See the individual methods pages for
information on whether and how this argument is used.}

\item{discard}{a string containing a method for discarding units outside a
region of common support. When a propensity score is estimated or supplied
to \code{distance} as a vector, the options are \code{"none"},
\code{"treated"}, \code{"control"}, or \code{"both"}. For \code{"none"}, no
units are discarded for common support. Otherwise, units whose propensity
scores fall outside the corresponding region are discarded. Can also be a
\code{logical} vector where \code{TRUE} indicates the unit is to be
discarded. Default is \code{"none"} for no common support restriction. See
Details.}

\item{reestimate}{if \code{discard} is not \code{"none"} and propensity
scores are estimated, whether to re-estimate the propensity scores in the
remaining sample. Default is \code{FALSE} to use the propensity scores
estimated in the original sample.}

\item{s.weights}{an optional numeric vector of sampling weights to be
incorporated into propensity score models and balance statistics. Can also
be specified as a string containing the name of variable in \code{data} to
be used or a one-sided formula with the variable on the right-hand side
(e.g., \code{~ SW}). Not all propensity score models accept sampling
weights; see \code{\link{distance}} for information on which do and do not,
and see \code{vignette("sampling-weights")} for details on how to use
sampling weights in a matching analysis.}

\item{replace}{for methods that allow it, whether matching should be done
with replacement (\code{TRUE}), where control units are allowed to be
matched to several treated units, or without replacement (\code{FALSE}),
where control units can only be matched to one treated unit each. See the
individual methods pages for information on whether and how this argument is
used. Default is \code{FALSE} for matching without replacement.}

\item{m.order}{for methods that allow it, the order that the matching takes
place. Allowable options depend on the matching method. The default of
\code{NULL} corresponds to \code{"largest"} when a propensity score is
estimated or supplied as a vector and \code{"data"} otherwise.}

\item{caliper}{for methods that allow it, the width(s) of the caliper(s) to
use in matching. Should be a numeric vector with each value named according
to the variable to which the caliper applies. To apply to the distance
measure, the value should be unnamed. See the individual methods pages for
information on whether and how this argument is used. Positive values require the distance between paired units to be no larger than the supplied caliper; negative values require the distance between paired units to be larger than the absolute value value of the supplied caliper. The default is \code{NULL} for no caliper.}

\item{std.caliper}{\code{logical}; when a caliper is specified, whether the
the caliper is in standard deviation units (\code{TRUE}) or raw units
(\code{FALSE}). Can either be of length 1, applying to all calipers, or of
length equal to the length of \code{caliper}. Default is \code{TRUE}.}

\item{ratio}{for methods that allow it, how many control units should be
matched to each treated unit in k:1 matching. Should be a single integer
value. See the individual methods pages for information on whether and how
this argument is used. The default is 1 for 1:1 matching.}

\item{verbose}{\code{logical}; whether information about the matching
process should be printed to the console. What is printed depends on the
matching method. Default is \code{FALSE} for no printing other than
warnings.}

\item{include.obj}{\code{logical}; whether to include any objects created in
the matching process in the output, i.e., by the functions from other
packages \code{matchit()} calls. What is included depends on the matching
method. Default is \code{FALSE}.}

\item{normalize}{\code{logical}; whether to rescale the nonzero weights in each treatment group to have an average of 1. Default is \code{TRUE}. See "How Matching Weights Are Computed" below for more details.}

\item{\dots}{additional arguments passed to the functions used in the
matching process. See the individual methods pages for information on what
additional arguments are allowed for each method.}
}
\value{
When \code{method} is something other than \code{"subclass"}, a
\code{matchit} object with the following components:

\item{match.matrix}{a matrix containing the matches. The row names correspond
to the treated units and the values in each row are the names (or indices)
of the control units matched to each treated unit. When treated units are
matched to different numbers of control units (e.g., with variable ratio matching or
matching with a caliper), empty spaces will be filled with \code{NA}. Not
included when \code{method} is \code{"full"}, \code{"cem"} (unless \code{k2k = TRUE}), \code{"exact"}, \code{"quick"}, or \code{"cardinality"} (unless \code{mahvars} is supplied and \code{ratio} is an integer).}
\item{subclass}{a factor
containing matching pair/stratum membership for each unit. Unmatched units
will have a value of \code{NA}. Not included when \code{replace = TRUE} or when \code{method = "cardinality"} unless \code{mahvars} is supplied and \code{ratio} is an integer.}
\item{weights}{a numeric vector of estimated matching weights. Unmatched and
discarded units will have a weight of zero.}
\item{model}{the fit object of
the model used to estimate propensity scores when \code{distance} is
specified as a method of estimating propensity scores. When
\code{reestimate = TRUE}, this is the model estimated after discarding
units.}
\item{X}{a data frame of covariates mentioned in \code{formula}, \code{exact}, \code{mahvars}, \code{caliper}, and \code{antiexact}.}
\item{call}{the \code{matchit()} call.}
\item{info}{information on the matching method and distance measures used.}
\item{estimand}{the argument supplied to \code{estimand}.}
\item{formula}{the \code{formula} supplied.}
\item{treat}{a vector of treatment status converted to zeros (0) and ones
(1) if not already in that format.}
\item{distance}{a vector of distance
values (i.e., propensity scores) when \code{distance} is supplied as a
method of estimating propensity scores or a numeric vector.}
\item{discarded}{a logical vector denoting whether each observation was
discarded (\code{TRUE}) or not (\code{FALSE}) by the argument to \code{discard}.}
\item{s.weights}{the vector of sampling weights supplied to the \code{s.weights} argument, if any.}
\item{exact}{a one-sided formula containing the variables, if any, supplied to \code{exact}.}
\item{mahvars}{a one-sided formula containing the variables, if any, supplied to \code{mahvars}.}
\item{obj}{when \code{include.obj = TRUE}, an object containing the intermediate results of the matching procedure. See
the individual methods pages for what this component will contain.}

When \code{method = "subclass"}, a \code{matchit.subclass} object with the same
components as above except that \code{match.matrix} is excluded and one
additional component, \code{q.cut}, is included, containing a vector of the
distance measure cutpoints used to define the subclasses. See
\code{\link{method_subclass}} for details.
}
\description{
\code{matchit()} is the main function of \emph{MatchIt} and performs
pairing, subset selection, and subclassification with the aim of creating
treatment and control groups balanced on included covariates. \emph{MatchIt}
implements the suggestions of Ho, Imai, King, and Stuart (2007) for
improving parametric statistical models by preprocessing data with
nonparametric matching methods.

This page documents the overall use of \code{matchit()}, but for specifics
of how \code{matchit()} works with individual matching methods, see the
individual pages linked in the Details section below.
}
\details{
Details for the various matching methods can be found at the following help
pages:
\itemize{
\item \code{\link{method_nearest}} for nearest neighbor matching
\item \code{\link{method_optimal}} for optimal pair matching
\item \code{\link{method_full}} for optimal full matching
\item \code{\link{method_quick}} for generalized (quick) full matching
\item \code{\link{method_genetic}} for genetic matching
\item \code{\link{method_cem}} for coarsened exact matching
\item \code{\link{method_exact}} for exact matching
\item \code{\link{method_cardinality}} for cardinality and profile matching
\item \code{\link{method_subclass}} for subclassification
}

The pages contain information on what the method does, which of the arguments above are
allowed with them and how they are interpreted, and what additional
arguments can be supplied to further tune the method. Note that the default
method with no arguments supplied other than \code{formula} and \code{data}
is 1:1 nearest neighbor matching without replacement on a propensity score
estimated using a logistic regression of the treatment on the covariates.
This is not the same default offered by other matching programs, such as
those in \emph{Matching}, \code{teffects} in Stata, or \verb{PROC PSMATCH}
in SAS, so care should be taken if trying to replicate the results of those
programs.

When \code{method = NULL}, no matching will occur, but any propensity score
estimation and common support restriction will. This can be a simple way to
estimate the propensity score for use in future matching specifications
without having to re-estimate it each time. The \code{matchit()} output with
no matching can be supplied to \code{summary()} to examine balance prior to
matching on any of the included covariates and on the propensity score if
specified. All arguments other than \code{distance}, \code{discard}, and
\code{reestimate} will be ignored.

See \code{\link{distance}} for details on the several ways to
specify the \code{distance}, \code{link}, and \code{distance.options}
arguments to estimate propensity scores and create distance measures.

When the treatment variable is not a \code{0/1} variable, it will be coerced
to one and returned as such in the \code{matchit()} output (see section
Value, below). The following rules are used: 1) if \code{0} is one of the
values, it will be considered the control and the other value the treated;
2) otherwise, if the variable is a factor, \code{levels(treat)[1]} will be
considered control and the other value the treated; 3) otherwise,
\code{sort(unique(treat))[1]} will be considered control and the other value
the treated. It is safest to ensure the treatment variable is a \code{0/1}
variable.

The \code{discard} option implements a common support restriction. It can
only be used when a distance measure is an estimated propensity score or supplied as a vector and is ignored for some matching
methods. When specified as \code{"treated"}, treated units whose distance
measure is outside the range of distance measures of the control units will
be discarded. When specified as \code{"control"}, control units whose
distance measure is outside the range of distance measures of the treated
units will be discarded. When specified as \code{"both"}, treated and
control units whose distance measure is outside the intersection of the
range of distance measures of the treated units and the range of distance
measures of the control units will be discarded. When \code{reestimate = TRUE} and \code{distance} corresponds to a propensity score-estimating
function, the propensity scores are re-estimated in the remaining units
prior to being used for matching or calipers.

Caution should be used when interpreting effects estimated with various
values of \code{estimand}. Setting \code{estimand = "ATT"} doesn't
necessarily mean the average treatment effect in the treated is being
estimated; it just means that for matching methods, treated units will be
untouched and given weights of 1 and control units will be matched to them
(and the opposite for \code{estimand = "ATC"}). If a caliper is supplied or
treated units are removed for common support or some other reason (e.g.,
lacking matches when using exact matching), the actual estimand targeted is
not the ATT but the treatment effect in the matched sample. The argument to
\code{estimand} simply triggers which units are matched to which, and for
stratification-based methods (exact matching, CEM, full matching, and
subclassification), determines the formula used to compute the
stratification weights.
\subsection{How Matching Weights Are Computed}{

Matching weights are computed in one of two ways depending on whether matching was done with replacement
or not.
\subsection{Matching without replacement and subclassification}{

For matching \emph{without} replacement (except for cardinality matching), including subclassification, each
unit is assigned to a subclass, which represents the pair they are a part of
(in the case of k:1 matching) or the stratum they belong to (in the case of
exact matching, coarsened exact matching, full matching, or
subclassification). The formula for computing the weights depends on the
argument supplied to \code{estimand}. A new "stratum propensity score"
(\eqn{p^s_i}) is computed for each unit \eqn{i} as \eqn{p^s_i = \frac{1}{n_s}\sum_{j: s_j =s_i}{I(A_j=1)}} where \eqn{n_s} is the size of subclass \eqn{s} and \eqn{I(A_j=1)} is 1 if unit \eqn{j} is treated and 0 otherwise. That is, the stratum propensity score for stratum \eqn{s} is the proportion of units in stratum \eqn{s} that are
in the treated group, and all units in stratum \eqn{s} are assigned that
stratum propensity score. This is distinct from the propensity score used for matching, if any. Weights are then computed using the standard formulas for
inverse probability weights with the stratum propensity score inserted:
\itemize{
\item for the ATT, weights are 1 for the treated
units and \eqn{\frac{p^s}{1-p^s}} for the control units
\item for the ATC, weights are
\eqn{\frac{1-p^s}{p^s}} for the treated units and 1 for the control units
\item for the ATE, weights are \eqn{\frac{1}{p^s}} for the treated units and \eqn{\frac{1}{1-p^s}} for the
control units.
}

For cardinality matching, all matched units receive a weight
of 1.
}

\subsection{Matching with replacement}{

For matching \emph{with} replacement, units are not assigned to unique strata. For
the ATT, each treated unit gets a weight of 1. Each control unit is weighted
as the sum of the inverse of the number of control units matched to the same
treated unit across its matches. For example, if a control unit was matched
to a treated unit that had two other control units matched to it, and that
same control was matched to a treated unit that had one other control unit
matched to it, the control unit in question would get a weight of \eqn{1/3 + 1/2 = 5/6}. For the ATC, the same is true with the treated and control labels
switched. The weights are computed using the \code{match.matrix} component
of the \code{matchit()} output object.
}

\subsection{Normalized weights}{

When \code{normalize = TRUE} (the default), in each treatment group, weights are divided by the mean of the nonzero
weights in that treatment group to make the weights sum to the number of
units in that treatment group (i.e., to have an average of 1).
}

\subsection{Sampling weights}{

If sampling weights are included through the
\code{s.weights} argument, they will be included in the \code{matchit()}
output object but not incorporated into the matching weights.
\code{\link[=match_data]{match_data()}}, which extracts the matched set from a \code{matchit} object,
combines the matching weights and sampling weights.
}

}
}
\examples{
data("lalonde")

# Default: 1:1 NN PS matching w/o replacement

m.out1 <- matchit(treat ~ age + educ + race + nodegree +
                    married + re74 + re75,
                  data = lalonde)
m.out1
summary(m.out1)

# 1:1 NN Mahalanobis distance matching w/ replacement and
# exact matching on married and race

m.out2 <- matchit(treat ~ age + educ + race + nodegree +
                    married + re74 + re75,
                  data = lalonde,
                  distance = "mahalanobis",
                  replace = TRUE,
                  exact = ~ married + race)
m.out2
summary(m.out2, un = TRUE)

# 2:1 NN Mahalanobis distance matching within caliper defined
# by a probit pregression PS

m.out3 <- matchit(treat ~ age + educ + race + nodegree +
                    married + re74 + re75,
                  data = lalonde,
                  distance = "glm",
                  link = "probit",
                  mahvars = ~ age + educ + re74 + re75,
                  caliper = .1,
                  ratio = 2)
m.out3
summary(m.out3, un = TRUE)

# Optimal full PS matching for the ATE within calipers on
# PS, age, and educ
\dontshow{if (requireNamespace("optmatch", quietly = TRUE)) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}
m.out4 <- matchit(treat ~ age + educ + race + nodegree +
                    married + re74 + re75,
                  data = lalonde,
                  method = "full",
                  estimand = "ATE",
                  caliper = c(.1, age = 2, educ = 1),
                  std.caliper = c(TRUE, FALSE, FALSE))
m.out4
summary(m.out4, un = TRUE)
\dontshow{\}) # examplesIf}
# Subclassification on a logistic PS with 10 subclasses after
# discarding controls outside common support of PS

s.out1 <- matchit(treat ~ age + educ + race + nodegree +
                    married + re74 + re75,
                  data = lalonde,
                  method = "subclass",
                  distance = "glm",
                  discard = "control",
                  subclass = 10)
s.out1
summary(s.out1, un = TRUE)
}
\references{
Ho, D. E., Imai, K., King, G., & Stuart, E. A. (2007). Matching
as Nonparametric Preprocessing for Reducing Model Dependence in Parametric
Causal Inference. \emph{Political Analysis}, 15(3), 199–236. \doi{10.1093/pan/mpl013}

Ho, D. E., Imai, K., King, G., & Stuart, E. A. (2011). MatchIt:
Nonparametric Preprocessing for Parametric Causal Inference. \emph{Journal of Statistical Software}, 42(8). \doi{10.18637/jss.v042.i08}
}
\seealso{
\code{\link[=summary.matchit]{summary.matchit()}} for balance assessment after matching, \code{\link[=plot.matchit]{plot.matchit()}} for plots of covariate balance and propensity score overlap after matching.
\itemize{
\item \code{vignette("MatchIt")} for an introduction to matching with \emph{MatchIt}
\item \code{vignette("matching-methods")} for descriptions of the variety of matching methods and options available
\item \code{vignette("assessing-balance")} for information on assessing the quality of a matching specification
\item \code{vignette("estimating-effects")} for instructions on how to estimate treatment effects after matching
\item \code{vignette("sampling-weights")} for a guide to using \emph{MatchIt} with sampling weights.
}
}
\author{
Daniel Ho, Kosuke Imai, Gary King, and Elizabeth Stuart wrote the original package. Starting with version 4.0.0, Noah Greifer is the primary maintainer and developer.
}
