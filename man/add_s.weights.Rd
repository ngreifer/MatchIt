% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/add_s.weights.R
\name{add_s.weights}
\alias{add_s.weights}
\title{Add sampling weights to a \code{matchit} object}
\usage{
add_s.weights(m, s.weights = NULL, data = NULL)
}
\arguments{
\item{m}{a \code{matchit} object; the output of a call to \code{\link[=matchit]{matchit()}},
typically with the \code{s.weights} argument unspecified.}

\item{s.weights}{an numeric vector of sampling weights to be added to the
\code{matchit} object. Can also be specified as a string containing the name
of variable in \code{data} to be used or a one-sided formula with the
variable on the right-hand side (e.g., \code{~ SW}).}

\item{data}{a data frame containing the sampling weights if given as a
string or formula. If unspecified, \code{add_s.weights()} will attempt to find
the dataset using the environment of the \code{matchit} object.}
}
\value{
a \code{matchit} object with an \code{s.weights} component
containing the supplied sampling weights. If \code{s.weights = NULL}, the original
\code{matchit} object is returned.
}
\description{
Adds sampling weights to a \code{matchit} object so that they are
incorporated into balance assessment and creation of the weights. This would
typically only be used when an argument to \code{s.weights} was not supplied
to \code{\link[=matchit]{matchit()}} (i.e., because they were not to be included in the estimation
of the propensity score) but sampling weights are required for generalizing
an effect to the correct population. Without adding sampling weights to the
\code{matchit} object, balance assessment tools (i.e., \code{\link[=summary.matchit]{summary.matchit()}}
and \code{\link[=plot.matchit]{plot.matchit()}}) will not calculate balance statistics correctly, and
the weights produced by \code{\link[=match_data]{match_data()}} and \code{\link[=get_matches]{get_matches()}} will not
incorporate the sampling weights.
}
\examples{

data("lalonde")

# Generate random sampling weights, just
# for this example
sw <- rchisq(nrow(lalonde), 2)

# NN PS match using logistic regression PS that doesn't
# include sampling weights
m.out <- matchit(treat ~ age + educ + race + nodegree +
                   married  + re74 + re75,
                 data = lalonde)

m.out

# Add s.weights to the matchit object
m.out <- add_s.weights(m.out, sw)

m.out #note additional output

# Check balance; note that sample sizes incorporate
# s.weights
summary(m.out, improvement = FALSE)

}
\seealso{
\code{\link[=matchit]{matchit()}}; \code{\link[=match_data]{match_data()}}
}
\author{
Noah Greifer
}
