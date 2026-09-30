\name{rerar_randomize}
\alias{rerar_randomize}
\title{Randomize the next enrollment batch}
\usage{
rerar_randomize(n, update, seed = NULL)
}
\arguments{
\item{n}{Number of participants in the next batch.}

\item{update}{An object returned by \code{\link[=rerar_update]{rerar_update()}}.}

\item{seed}{Optional random seed.}
}
\value{
A data frame containing binary assignments and their common stored
assignment probability.
}
\description{
Randomize the next enrollment batch
}
