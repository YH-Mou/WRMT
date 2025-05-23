% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/WRMT.R
\name{WRMT.comparison}
\alias{WRMT.comparison}
\title{Pairwise comparisons for WRAT}
\usage{
WRMT.comparison(
  data,
  threshold,
  stage = NULL,
  direction = NULL,
  strata = NULL,
  strata.weight = NULL
)
}
\arguments{
\item{data}{A data frame with each row representing a participant.
The columns should be orders as:
\itemize{
\item \code{ID}: Participant ID, must be numeric and unique.
\item \code{Treatment Group Indicator}: 1 (True) = Treatment, 0 (False) = Control.
\item \code{Outcome 1}: Observed outcome 1 in numeric format.
\item \code{Censor Status of Outcome 1}: 1 (True) = Censored, 0 (False) = Not Censored.
\item \code{Outcome 2}: Observed outcome 2 in numeric format.
\item \code{Censor Status of Outcome 2}: 1 (True) = Censored, 0 (False) = Not Censored.
\item similar for outcome 3, 4, and more.
}}

\item{threshold}{A numeric vector indicating thresholds for each stage.}

\item{stage}{A numeric vector indicating the outcome measure considered at
each stage of comparison (1 for outcome 1, 2 for outcome 2, etc.).
If specified, the length must be the same with the length of threshold parameter.
Default to alternating across outcomes when the length of threshold is a multiply of
number of outcomes. E.g., with two outcomes and four thresholds specified,
the default stage is (1,2,1,2).}

\item{direction}{A numeric vector the more favorable results of each outcome.
1 for favoring larger outcome values (e.g., surviving longer time) and -1 for favoring
smaller outcome values (e.g., having fewer adverse events).
efault to favoring larger values for all outcomes.}

\item{strata}{A numeric vector indicates the stratum that each participant belongs to.
Only needed for analysis with stratification. Default to leaving blank for non-stratified
analysis.}

\item{strata.weight}{A data frame indicates weights of each stratum.
Default to weighting base on numbers of participants within strata.
If specified, the columns should be orders as:
\itemize{
\item \code{stratum}: Stratum index, should be numeric,
unique and corresponds to \code{strata}.
\item \code{weight}: Weight for each stratum, should be numeric.
}}
}
\value{
A list of comparison results:
\describe{
\item{pair.data.list}{A list that contains data frames of pairwise comparison
results and data used for comparison. Each stratum corresponds to one data frame.}
\item{stage.results.list}{A list that contains matrices of comparison results of each
stage. Each stratum corresponds to one matrix.
1 for win, -1 for loss, 0 for tie, and NA for results determined by previous stages.}
\item{treatment}{A numeric vector indicates the treatment status.}
\item{strata}{A numeric vector indicates the stratum that each participant belongs to.}
\item{strata.weight}{A data frame contains weighting information
for stratification analysis.}
}
}
\description{
This function performs pairwise comparisons for WRMT analysis. To be called
by function \code{WRMT()}.
}
\examples{
# Example 1: stratified comparison
strata.compare <- WRMT.comparison(data = sample_data[,1:6],
                                  threshold = c(0,0), stage = c(1,2),
                                  direction = c(1,1),
                                  strata = sample_data$strata,
                                  strata.weight = NULL)

# Example 2: non-stratified comparison
nonstrata.compare <- WRMT.comparison(data = sample_data[,1:6],
                                     threshold = c(0,0), stage = c(1,2),
                                     direction = c(1,1),
                                     strata.weight = NULL)


}
