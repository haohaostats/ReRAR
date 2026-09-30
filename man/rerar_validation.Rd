\name{rerar_validation}
\alias{rerar_validation}
\title{Create external paired-validation counts}
\usage{
rerar_validation(arm, se_success, se_failure, sp_success, sp_failure)
}
\arguments{
\item{arm}{Treatment indicator. Supply one row for each of arms 0 and 1.}

\item{se_success, se_failure}{Counts among final-outcome successes for which
the provisional outcome is respectively positive or negative.}

\item{sp_success, sp_failure}{Counts among final-outcome failures for which
the provisional outcome is respectively negative or positive.}
}
\value{
A data frame of class \code{rerar_validation}.
}
\description{
Creates the treatment-specific sensitivity and specificity counts used by
the revision likelihood. Arm 0 is the control treatment and arm 1 is the
experimental treatment.
}
