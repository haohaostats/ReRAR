\name{exam_example}
\alias{exam_example}
\docType{data}
\title{EXAM-based worked example of revision-aware randomization}
\description{
One case-study object combining published paired-assessment counts with a
synthetic trial generated for demonstration. Individual participant records,
enrollment and outcome-return times are simulated, not original EXAM records.
}
\usage{data(exam_example)}
\format{
A named list with five components:
\describe{
\item{validation}{Two-row \code{rerar_validation} object, with arm 0 denoting
placebo and arm 1 cabozantinib. The four count columns are sensitivity
successes/failures and specificity successes/failures.}
\item{trial}{A \code{rerar_history} object for 154 participants, frozen before
batch 8. Columns \code{id}, \code{A}, \code{Q}, \code{QZ}, \code{R},
\code{RY} encode participant ID, treatment, provisional visibility and outcome,
and final visibility and outcome. Unavailable outcomes are \code{NA}.}
\item{final}{330 synthetic participant records. Columns \code{id},
\code{entry_month}, \code{A}, \code{Z}, \code{Y}, \code{pi},
\code{Z_available}, \code{Y_available}, \code{working0}, \code{working1}
contain identifiers, enrollment batch, treatment, complete outcomes, stored
assignment probability, outcome-return months and prespecified working means.}
\item{trajectory}{15 batch updates: \code{month}, \code{probability},
\code{target}, \code{cumulative_allocation}, \code{provisional_observed},
\code{final_observed}. Outcome counts refer to the freeze before each batch;
cumulative allocation refers to the end of the batch.}
\item{metadata}{Provenance, seed, trial settings and the month-8 snapshot definition.}
}}
\details{
The synthetic trial uses final success probabilities 55/111 and 150/219,
with 15 batches of 22 participants. Provisional outcomes return two months
after enrollment; an additional log-normal delay has median 3.05 months and
90th percentile 7.5 months. Allocation bounds are 0.20 and 0.78 and DBCD gamma
is 2. Both final-analysis working means are fixed at 0.5 before enrollment.
The seed is 20260930. This single realization does not establish operating
characteristics or reproduce the manuscript's 10,000-replicate results.
}
\source{FDA statistical review for cabozantinib, Table 13:
\url{https://www.accessdata.fda.gov/drugsatfda_docs/nda/2012/203756Orig1s000StatR.pdf}}
\examples{
data(exam_example)
x <- exam_example
update <- rerar_update(x$trial, x$validation, lower = 0.20, upper = 0.78)
print(update)
with(x$final, rerar_analyze(A, Y, pi, working0, working1))
}
\keyword{datasets}
