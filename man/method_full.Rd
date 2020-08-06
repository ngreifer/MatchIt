\name{method_full}
\alias{method_full}

\title{Optimal Full Matching}

\description{
In \code{\link{matchit}}, setting \code{method = "full"} performs optimal full matching, which is a form of subclassification wherein all units, both treatment and control (i.e., the "full" sample), are assigned to a subclass and receive at least one match. The matching is optimal in the sense that that sum of the absolute distances between the treated and control units in each subclass are as small as possible. The method relies on and is a wrapper for the \code{\link[optmatch]{fullmatch}} function in the \pkg{optmatch} package. Advantages of optimal full matching include that the matching order is not required to be specified, units do not need to be discarded, and it is less likely that extreme within-subclass distances will be large, unlike with standard subclassification. The primary output of full matching is a set of matching weights that can be applied to the matched sample; in the way, full matching can be seen as a robust alternative to propensity score weighting, robust in the sense that the propensity score model does not need to be correct for good matches to be found.

This page details the allowable arguments with \code{method = "fullmatch"}. See \code{\link{matchit}} for an explanation of what each argument means in a general context and how it can be specified.
}
\usage{
matchit(formula, data = NULL, method = "full", distance = "glm",
        link = "logit", distance.options=list(), ratio = 1,
        caliper = NULL, exact = NULL, mahvars = NULL,
        discard = "none", reestimate = FALSE, verbose = FALSE,
        ...)

}

\arguments{
  \item{formula}{
a two-sided \code{\link{formula}} object containing the treatment and covariates to be used in creating the distance measure used in the matching. This formula will be supplied to the functions that estimate the distance measure.
}
  \item{data}{
a data frame containing the variables named in \code{formula}. If not found in \code{data}, the variables will be saught in the environment.
}
  \item{method}{
set here to \code{"full"}.
}
  \item{distance}{
the distance measure to be used. See \code{\link{distance}} for allowable options. When set to \code{"mahalanobis"}, optimal full Mahalanobis distance matching will be performed on the variables named in \code{formula}.
}
  \item{link}{
when \code{distance} is specified as a string and not \code{"mahalanobis"}, an additional argument controlling the link function used in estimating the distance measure. See \code{\link{distance}} for allowable options with each option.
}
  \item{distance.options}{
a named list containing additional arguments supplied to the function that estimates the distance measure as determined by the argument to \code{distance}.
}
  \item{caliper}{
the width of the caliper used for caliper matching in units of the standard deviation of the distance measure. Calipers are processed using \code{\link[optmatch]{caliper}} in \pkg{optmatch}. Not allowed when \code{distance = "mahalanobis"}. See Notes.
}
  \item{exact}{
for which variables exact matching should take place. Exact matching is processed using \code{\link[optmatch]{exactMatch}} in \pkg{optmatch}.
}
  \item{mahvars}{
for which variables Mahalanobis matching should take place when a distance measure other than \code{"mahalanobis"} is used (e.g., for caliper matching or to discard units for common support). If specified, the distance measure will not be used in matching.
}
  \item{discard}{
a string containing a method for discarding units outside a region of common support. Only allowed when \code{distance} is not \code{"mahalanobis"}.
}
  \item{reestimate}{
if \code{discard} is not \code{"none"}, whether to re-estimate the propensity score in the remaining sample prior to matching.
}
  \item{verbose}{
\code{logical}; whether information about the matching process should be printed to the console.
}
  \item{\dots}{
additional arguments passed to \code{\link[optmatch]{fullmatch}} in \pkg{optmatch}. Potentially useful arguments include \code{min.controls}, \code{max.controls}, \code{omit.fraction}, \code{mean.controls}, and \code{tol}. See the \code{\link[optmatch]{fullmatch}} documentation for details.
}

The arguments \code{replace}, \code{ratio}, and \code{m.order} are ignored with a warning.
}
\notes{
Due to a bug in \pkg{optmatch} (version 0.9-13), a caliper can only be used when \code{min.controls}, \code{max.controls}, \code{omit.fraction}, and \code{mean.controls} are left at their defaults.
}
\section{Outputs}{
All outputs described in \code{\link{matchit}} are returned with \code{method = "full"} except for \code{match.matrix}. This is because matching strata are not indexed by treated units as they are in some other forms of matching.
}
\seealso{
\code{\link{matchit}} for a detailed explanation of the inputs and outputs of a call to \code{matchit}.

\code{\link[optmatch]{fullmatch}} in \pkg{optmatch}, which is the workhorse.

\code{\link{method_optimal}} for optimal pair matching, which is a special case of optimal full matching, and which relies on similar machinery. Results from \code{method = "optimal"} can be replicated with \code{method = "full"} by setting \code{min.controls}, \code{max.controls}, and \code{mean.controls} to \code{ratio}.
}
\details{
Mahalanobis distance matching can be done one of two ways:

1) If a propensity score needs to be estimated for any reason, e.g., for common support with \code{discard} or for creating a caliper, \code{distance} should be whatever method is used to estimate the propensity score or a vector of distance measures, i.e., it should not be \code{"mahalanobis"}. Use \code{mahvars} to specify the variables used to create the Mahalanobis distance. For example, to perform Mahalanobis within a propensity score caliper, the following could be run:

\code{matchit(treat ~ X1 + X2 + X3, method = "full", distance = "glm",
              caliper = .25, mahvars = ~ X1 + X2)}

With this code, \code{X1}, \code{X2}, and \code{X3} are used to estimate the propensity score (using the \code{"glm"} method, which by default is logistic regression), which is used to create a matching caliper. The actual matching occurs on the Mahalanobis distance computed only using \code{X1} and \code{X2}, which are supplied to \code{mahvars}. Units whose propensity score difference is larger than the caliper will not be paired, and some treated units may therefore not receive a match. The estimated propensity scores will be included in the \code{distance} component of the \code{matchit} output.

2) If no propensity score needs to be estimated, \code{distance} should be set to \code{"mahalanobis"}, and Mahalanobis distance matching will occur on all the variables in \code{formula}. Arguments to \code{caliper}, \code{discard}, and \code{mahvars} will be ignored. For example, to perform simple Mahalanobis distance matching, the following could be run:

\code{matchit(treat ~ X1 + X2, method = "full", distance = "mahalanobis")}

With this code, the Mahalanobis distance is computed using \code{X1} and \code{X2}, and matching occurs on this distance. The \code{distance} component of the \code{matchit} output will be empty.

}
\references{
In a manuscript, be sure to cite the following paper if using \code{matchit} with \code{method = "full"}:

Hansen, B. B., & Klopfer, S. O. (2006). Optimal Full Matching and Related Designs via Network Flows. Journal of Computational and Graphical Statistics, 15(3), 609–627. \doi{10.1198/106186006X137047}

For example, a sentence might read:

\emph{Optimal full matching was performed using the MatchIt package (Ho, Imai, King, & Stuart, 2011) in R, which calls functions from the optmatch package (Hansen & Klopfer, 2006).}

Theory is also developed in the following article:

Hansen, B. B. (2004). Full Matching in an Observational Study of Coaching for the SAT. Journal of the American Statistical Association, 99(467), 609–618. \doi{10.1198/016214504000000647}

}
\examples{
data("lalonde")

# Optimal full PS matching
m.out1 <- matchit(treat ~ age + educ + race + nodegree +
                    married + re74 + re75, data = lalonde,
                  method = "full")
summary(m.out1)

# Optimal full Mahalanobis distance matching within a PS caliper
m.out2 <- matchit(treat ~ age + educ + race + nodegree +
                    married + re74 + re75, data = lalonde,
                  method = "full", caliper = .01,
                  mayvars = ~ age + educ + re74 + re75)
summary(m.out2)
}