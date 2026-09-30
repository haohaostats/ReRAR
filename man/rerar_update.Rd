\name{rerar_update}
\alias{rerar_update}
\title{Compute the next ReRAR allocation probability}
\usage{
rerar_update(history, validation, gamma = 2, lower = 0.39, upper = 0.62)
}
\arguments{
\item{history}{A trial history created by \code{\link[=rerar_history]{rerar_history()}}.}

\item{validation}{External paired counts created by \code{\link[=rerar_validation]{rerar_validation()}}.}

\item{gamma}{DBCD tracking parameter.}

\item{lower, upper}{Prespecified allocation boundaries.}
}
\value{
An object of class \code{rerar_update} with the next-batch probability,
target allocation, evidence statistic, and estimated final response rates.
}
\description{
Compute the next ReRAR allocation probability
}
