\name{rerar_analyze}
\alias{rerar_analyze}
\title{Analyze the final binary outcome from a ReRAR trial}
\usage{
rerar_analyze(
  A,
  Y,
  pi,
  working0 = 0.5,
  working1 = 0.5,
  alternative = c("greater", "less", "two.sided"),
  conf.level = 0.95
)
}
\arguments{
\item{A}{Binary treatment assignment.}

\item{Y}{Observed final binary outcome.}

\item{pi}{Stored probability of assignment to the experimental arm.}

\item{working0, working1}{Fixed or history-predictable working means. Scalar
values are expanded to all participants. The default is 0.5 in both arms.}

\item{alternative}{One of \code{"greater"}, \code{"less"}, or \code{"two.sided"}.}

\item{conf.level}{Confidence level for the Wald interval.}
}
\value{
An object of class \code{rerar_analysis}.
}
\description{
Uses the stored assignment probability for each randomized participant and
the contrast-stabilized augmented inverse-probability weighted estimator.
}
