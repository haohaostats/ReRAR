\name{rerar_history}
\alias{rerar_history}
\title{Create an observed ReRAR trial history}
\usage{
rerar_history(
  A,
  Z = rep(NA_integer_, length(A)),
  Y = rep(NA_integer_, length(A)),
  id = seq_along(A)
)
}
\arguments{
\item{A}{Binary treatment assignment, where 1 is experimental treatment.}

\item{Z}{Provisional binary outcome. Use \code{NA} when it is not yet visible.}

\item{Y}{Final binary outcome. Use \code{NA} when it has not yet returned.}

\item{id}{Optional participant identifier.}
}
\value{
A data frame of class \code{rerar_history} containing the internal
observability indicators required by \code{\link[=rerar_update]{rerar_update()}}.
}
\description{
Create an observed ReRAR trial history
}
