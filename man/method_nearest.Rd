\name{method_nearest}
\alias{method_nearest}
\title{
Nearest Neighbor Matching
}
\description{
In \code{\link{matchit}}, setting \code{method = "nearest"} performs greedy nearest neighbor matching. A distance is computed between each treated unit and each control unit, and, one by one, each treated unit is assigned a control unit as a match. The matching is "greedy" in the snese that there is no action taken to optimize an overall criterion; each match is selected without considering the other matches that may occur subsequently. Nearest neighbor matching implemented here does not rely on any other package.

This page details the allowable arguments with \code{method = "nearest"}. See \code{\link{matchit}} for an explanation of what each argument means in a general context and how it can be specified.
}
\usage{
matchit(formula, data = NULL, method = "genetic", distance = "glm",
        link = "logit", distance.options=list(), replace = FALSE,
        ratio = 1, caliper = NULL, m.order = NULL, exact = NULL,
        mahvars = NULL, discard = "none", reestimate = FALSE,
        verbose = FALSE, ...)
}
\arguments{
   \item{formula}{
a two-sided \code{\link{formula}} object containing the treatment and covariates to be used in creating the distance measure used in the matching.
}
  \item{data}{
a data frame containing the variables named in \code{formula}. If not found in \code{data}, the variables will be saught in the environment.
}
  \item{method}{
set here to \code{"nearest"}.
}
  \item{distance}{
the distance measure to be used. See \code{\link{distance}} for allowable options.
}
  \item{link}{
when \code{distance} is specified as a string and not \code{"mahalanobis"}, an additional argument controlling the link function used in estimating the distance measure. See \code{\link{distance}} for allowable options with each option.
}
  \item{distance.options}{
a named list containing additional arguments supplied to the function that estimates the distance measure as determined by the argument to \code{distance}.
}
  \item{replace}{
whether matching should be done with replacement.
}
  \item{ratio}{
how many control units should be matched to each treated unit for k:1 matching. Should be a single integer value.
}
  \item{caliper}{
the width of the caliper used for caliper matching in units of the standard deviation of the distance measure. Not allowed when \code{distance = "mahalanobis"}.
}
  \item{m.order}{
the order that the matching takes place. The default for \code{distance = "mahalanobis"} is \code{"data"}. Otherwise, the default is \code{"largest"}. See \code{\link{matchit}} for allowable options.
}
  \item{exact}{
for which variables exact matching should take place.
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
  \item{estimand}{
a string containing the desired estimand. Allowable options include \code{"ATT"} and \code{"ATC"}. See Details.
}
  \item{verbose}{
\code{logical}; whether information about the matching process should be printed to the console.
}
  \item{\dots}{
ignored.
}
}
\section{Outputs}{
All outputs described in \code{\link{matchit}} are returned with \code{method = "nearest"}. When \code{repalce = TRUE}, the \code{subclass} component is omitted.
}
\details{
Mahalanobis distance matching can be done one of two ways:

1) If a propensity score needs to be estimated for any reason, e.g., for common support with \code{discard} or for creating a caliper, \code{distance} should be whatever method is used to estimate the propensity score or a vector of distance measures, i.e., it should not be \code{"mahalanobis"}. Use \code{mahvars} to specify the variables used to create the Mahalanobis distance. For example, to perform Mahalanobis within a propensity score caliper, the following could be run:

\code{matchit(treat ~ X1 + X2 + X3, method = "nearest", distance = "glm",
              caliper = .25, mahvars = ~ X1 + X2)}

With this code, \code{X1}, \code{X2}, and \code{X3} are used to estimate the propensity score (using the \code{"glm"} method, which by default is logistic regression), which is used to create a matching caliper. The actual matching occurs on the Mahalanobis distance computed only using \code{X1} and \code{X2}, which are supplied to \code{mahvars}. Units whose propensity score difference is larger than the caliper will not be paired, and some treated units may therefore not receive a match. The estimated propensity scores will be included in the \code{distance} component of the \code{matchit} output.

2) If no propensity score needs to be estimated, \code{distance} should be set to \code{"mahalanobis"}, and Mahalanobis distance matching will occur on all the variables in \code{formula}. Arguments to \code{caliper}, \code{discard}, and \code{mahvars} will be ignored. For example, to perform simple Mahalanobis distance matching, the following could be run:

\code{matchit(treat ~ X1 + X2, method = "nearest", distance = "mahalanobis")}

With this code, the Mahalanobis distance is computed using \code{X1} and \code{X2}, and matching occurs on this distance. The \code{distance} component of the \code{matchit} output will be empty.

\subsection{Estimand}{
The \code{estimand} argument controls whether control units are selected to be matched with treated units (\code{estimand = "ATT"}) or treated units are selected to be matched with control units (\code{estimand = "ATC"}). The "focal" group (e.g., the treated units for the ATT) is typically made to be the smaller treatmentgroup, and a warning will be thrown if it is not set that way unless \code{replace = TRUE}. Setting \code{estimand = "ATC"} is equivalent to swapping all treated and control labels for the treatment variable. When \code{estimand = "ATC"}, the default \code{m.order} is \code{"smallest"}, and the \code{match.matrix} component of the output will have the names of the control units as the rownames and be filled with the names of the matched treated units (opposite to when \code{estimand = "ATT"}). Note that the argument supplied to \code{estimand} doesn't necessarily correspond to the estimand actually targeted; it is merely a switch to trigger which treatment group is considered "focal".
}

}
\references{
In a manuscript, you don't need to cite another package when using \code{method = "nearest"} because the matching is performed completely within \pkg{MatchIt}. For example, a sentence might read:

\emph{Nearest neighbor matching was performed using the MatchIt package (Ho, Imai, King, & Stuart, 2011) in R.}
}

\seealso{
\code{\link{matchit}} for a detailed explanation of the inputs and outputs of a call to \code{matchit}.

\code{\link{method_optimal}} for optimal pair matching, which is similar to greedy nearest neighbor matching except that an overall distance criterion is minimized.
}
\examples{
data("lalonde")

# 1:1 greedy NN matching on the PS
m.out1 <- matchit(treat ~ age + educ + race + nodegree +
                    married + re74 + re75, data = lalonde,
                  method = "nearest")
summary(m.out1)

# 2:1 genetic matching with replacement without PS
m.out2 <- matchit(treat ~ age + educ + race + nodegree +
                    married + re74 + re75, data = lalonde,
                  method = "genetic", replace = TRUE,
                  ratio = 2, distance = "mahalanobis")
summary(m.out2)
}