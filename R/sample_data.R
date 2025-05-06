#' Sample dataset for WRMT analysis
#'
#' A simulated dataset for demonstration purposes.
#' This dataset has to time-to-event endpionts with poential administrative
#' right censoring (i.e., maximal follow-up time set to 1000 time units),
#' where the second event may be censored by the first event
#' (i.e., semi-competing risk).
#' The columns include:
#' \itemize{
#'   \item \code{ID}: Numeric participant ID.
#'   \item \code{treat}: 1 = Treatment group, 0 = Control group.
#'   \item \code{T1}: Observed value of the frist time-to-event endpoint.
#'   \item \code{C1}: Censoring indicator for T1 (1 = Censored, 0 = Event observed).
#'   \item \code{T2}: Observed value of the second time-to-event endpoint.
#'   \item \code{C2}: Censoring indicator for T2.
#'   \item \code{strata}: Strata information (with 3 different strata).
#' }
#'
#' @format A data frame with 2000 rows and 7 columns.
#' @examples
#' data(sample_data)
#' head(sample_data)
"sample_data"
