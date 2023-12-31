% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/validate_pairs.R
\name{validate_pairs}
\alias{validate_pairs}
\title{Validate sample pairing for pair-ended RNA seq}
\usage{
validate_pairs(pairs, rows, reps)
}
\arguments{
\item{pairs}{Numeric vector of integers and/or \code{NA}. Positions with
the same integer value indicate samples that are paired in a pair-ended
experiment.}

\item{rows}{Numeric vector of integers indicating selected rows.}

\item{reps}{Numeric vector of integers and/or \code{NA}. Positions with the
same integer value indicate samples that replicates.}
}
\value{
\code{TRUE} if the pairing is valid, otherwise \code{FALSE}.
}
\description{
Used internally by \link{select_pairs} and externally by drugseqr
}
\examples{
# example if valid pairing
pairs <- rep(NA, 4)
rows <- c(1, 2)
reps <- rep(NA, 4)
validate_pairs(pairs, rows, reps)
}
