#' WRMT Analysis
#' @description
#' Perform WRMT (win ratio with multiple thresholds) analysis.
#'
#' @param data A data frame with each row representing a participant.
#' The columns should be orders as:
#' \itemize{
#'   \item \code{ID}: Participant ID, must be numeric and unique.
#'   \item \code{Treatment Group Indicator}: 1 (True) = Treatment, 0 (False) = Control.
#'   \item \code{Outcome 1}: Observed outcome 1 in numeric format.
#'   \item \code{Censor Status of Outcome 1}: 1 (True) = Censored, 0 (False) = Not Censored.
#'   \item \code{Outcome 2}: Observed outcome 2 in numeric format.
#'   \item \code{Censor Status of Outcome 2}: 1 (True) = Censored, 0 (False) = Not Censored.
#'   \item similar for outcome 3, 4, and more.
#' }
#' @param threshold A numeric vector indicating thresholds for each stage.
#' @param stage A numeric vector indicating the outcome measure considered at
#' each stage of comparison (1 for outcome 1, 2 for outcome 2, etc.).
#' If specified, the length must be the same with the length of threshold parameter.
#' Default to alternating across outcomes when the length of threshold is a multiply of
#' number of outcomes. E.g., with two outcomes and four thresholds specified,
#' the default stage is (1,2,1,2).
#' @param direction A numeric vector the more favorable results of each outcome.
#' 1 for favoring larger outcome values (e.g., surviving longer time) and -1 for favoring
#' smaller outcome values (e.g., having fewer adverse events).
#' efault to favoring larger values for all outcomes.
#' @param strata A numeric vector indicates the stratum that each participant belongs to.
#' Only needed for analysis with stratification. Default to leaving blank for non-stratified
#' analysis.
#' @param strata.weight A data frame indicates weights of each stratum.
#' Default to weighting base on numbers of participants within strata.
#' If specified, the columns should be orders as:
#' \itemize{
#'   \item \code{stratum}: Stratum index, should be numeric,
#'   unique and corresponds to \code{strata}.
#'   \item \code{weight}: Weight for each stratum, should be numeric.
#' }
#'
#' @returns A list of analysis results:
#' \describe{
#'   \item{WR}{A numeric value. Estimated win ratio for treatment versus control.}
#'   \item{P}{A numeric value. Two-sided P value of the Finkelstein-Schoenfeld (FS) test.}
#'   \item{Stage.decompose}{A list of stage-level comparison results for each stratum.
#'   Each stratum corresponds to a matrix. A single matrix will be returned for
#'   non-stratified analysis for simplicity. The decomposition matrix presents
#'   proportions of wins, ties, and losses (obtained by the treatment group) at each stage.
#'   Proportions are always obtained with respect to the total number of comparisons between
#'   treatment and control groups.}
#'   \item{strata.summary}{A data frame presents strata weights and stratum-specific
#'   WR results.}
#'   \item{FS.stat}{A numeric value of the FS test statistic.}
#'   \item{FS.stat.var}{A numeric value of the variance of the FS test statistic.}
#' }
#' @export
#'
#' @examples
#' # Example 1: stratified WRMT analysis with two time-to-event endpoints.
#' #   Stages: endpoint 1 with threshold 10, endpoint 2 with threshold 7,
#' #   endpoint 1 with thresholds 0, and endpoint 2 with threshold 0.
#' #   That is stratified WRMT(10,7,0,0) analysis.
#'
#' analysis.results <- WRMT(data = sample_data[,1:6],
#'                          threshold = c(0,0), stage = c(1,2),
#'                          direction = c(1,1),
#'                          strata = NULL,
#'                          strata.weight = NULL)
#' cat("WR:",analysis.results$WR, "; P-value:", analysis.results$P,"\n")
#'
#'
#' # Example 2: non-stratified WRMT analysis with two time-to-event endpoints.
#' #   Stages: endpoint 1 with threshold 10, endpoint 2 with threshold 7,
#' #   endpoint 1 with thresholds 0, and endpoint 2 with threshold 0.
#' #   That is non-stratified WRMT(10,7,0,0) analysis.
#'
#' analysis.results <- WRMT(data = sample_data[,1:6],
#'                          threshold = c(0,0), stage = c(1,2),
#'                          direction = c(1,1),
#'                          strata = NULL,
#'                          strata.weight = NULL)
#' cat("WR:",analysis.results$WR, "; P-value:", analysis.results$P,"\n")
#'
#' # Example 3: stratified standard Win Ratio analysis with
#' #   two time-to-event endpoints
#' analysis.results <- WRMT(data = sample_data[,1:6],
#'                          threshold = c(0,0), stage = c(1,2),
#'                          direction = c(1,1),
#'                          strata = sample_data$strata,
#'                          strata.weight = NULL)
#' cat("WR:",analysis.results$WR, "; P-value:", analysis.results$P,"\n")
#'
#' # Example 4: non-stratified standard Win Ratio analysis with
#' #   two time-to-event endpoints
#' analysis.results <- WRMT(data = sample_data[,1:6],
#'                          threshold = c(0,0), stage = c(1,2),
#'                          direction = c(1,1),
#'                          strata = NULL,
#'                          strata.weight = NULL)
#' cat("WR:",analysis.results$WR, "; P-value:", analysis.results$P,"\n")
#'
WRMT <- function (data, threshold, stage = NULL, direction = NULL,
                  strata = NULL, strata.weight = NULL) {
  ## check input information

  ## perform analysis and output results
  intermediate.res <- WRMT.comparison(data = data, threshold = threshold,
                                      stage = stage, direction = direction,
                                      strata = strata, strata.weight = strata.weight)
  res.output <- comparison.result(intermediate.res)
  return(res.output)
}




#' WRAT Analysis
#'
#' @description
#' Perform WRAT (win ratio with adaptive thresholds) analysis.
#'
#'
#' @param data A data frame with each row representing a participant.
#' The columns should be orders as:
#' \itemize{
#'   \item \code{ID}: Participant ID, must be numeric and unique.
#'   \item \code{Treatment Group Indicator}: 1 (True) = Treatment, 0 (False) = Control.
#'   \item \code{Outcome 1}: Observed outcome 1 in numeric format.
#'   \item \code{Censor Status of Outcome 1}: 1 (True) = Censored, 0 (False) = Not Censored.
#'   \item \code{Outcome 2}: Observed outcome 2 in numeric format.
#'   \item \code{Censor Status of Outcome 2}: 1 (True) = Censored, 0 (False) = Not Censored.
#'   \item similar for outcome 3, 4, and more.
#' }
#' @param method Method used to obtain adaptive thresholds.
#' "weights" or "calipers" as detailed below.
#' @param weights A numeric vector indicates the weights when method "weights" is employed.
#' Default to 1 for all outcomes.
#' @param calipers A numeric value (or vector) indicates the calipers.
#'
#' For \code{method="weights"}, a single value that will be applied to all outcomes.
#' Default to 0.2.
#'
#' For \code{method="calipers"}, a vector that contains calipers for each outcome.
#' Input required for \code{method="calipers"}.
#'
#' @param clinical.threshold
#' A numeric vector represents the minimal clinical thresholds.
#' Default to 0 for all endpoints.
#'
#' @param direction A numeric vector the more favorable results of each outcome.
#' 1 for favoring larger outcome values (e.g., surviving longer time) and -1 for favoring
#' smaller outcome values (e.g., having fewer adverse events).
#' @param strata A numeric vector indicates the stratum that each participant belongs to.
#' Only needed for analysis with stratification. Default to leaving blank for non-stratified
#' analysis.
#' @param strata.weight A data frame indicates weights of each stratum.
#' Default to weighting base on numbers of participants within strata.
#' If specified, the columns should be orders as:
#' \itemize{
#'   \item \code{stratum}: Stratum index, should be numeric,
#'   unique and corresponds to \code{strata}.
#'   \item \code{weight}: Weight for each stratum, should be numeric.
#' }
#'
#' @returns A list of analysis results:
#' \describe{
#'   \item{WR}{A numeric value. Estimated win ratio for treatment versus control.}
#'   \item{P}{A numeric value. Two-sided P value of the Finkelstein-Schoenfeld (FS) test.}
#'   \item{Stage.decompose}{A list of stage-level comparison results for each stratum.
#'   Each stratum corresponds to a matrix. A single matrix will be returned for
#'   non-stratified analysis for simplicity. The decomposition matrix presents
#'   proportions of wins, ties, and losses (obtained by the treatment group) at each stage.
#'   Proportions are always obtained with respect to the total number of comparisons between
#'   treatment and control groups.}
#'   \item{strata.summary}{A data frame presents strata weights and stratum-specific
#'   WR results.}
#'   \item{FS.stat}{A numeric value of the FS test statistic.}
#'   \item{FS.stat.var}{A numeric value of the variance of the FS test statistic.}
#' }
#'
#' @export
#'
#' @examples
#' # Example 1: stratified WRAT analysis with two time-to-event endpoints.
#' #   The "weights" method employed.
#' analysis.results <- WRAT(data = sample_data[,1:6],
#'                          method = "weights",
#'                          weights = c(1,1), calipers = 0.2,
#'                          clinical.threshold = c(0,0),
#'                          direction = c(1,1),
#'                          strata = sample_data$strata,
#'                          strata.weight = NULL)
#' cat("WR:",analysis.results$WR, "; P-value:", analysis.results$P,"\n")
#'
#' # Example 2: stratified WRAT analysis with two time-to-event endpoints.
#' #   The "calipers" method employed.
#' analysis.results <- WRAT(data = sample_data[,1:6],
#'                          method = "calipers",
#'                          weights = NULL, calipers = c(0.2,0.3),
#'                          clinical.threshold = c(0,0),
#'                          direction = c(1,1),
#'                          strata = sample_data$strata,
#'                          strata.weight = NULL)
#' cat("WR:",analysis.results$WR, "; P-value:", analysis.results$P,"\n")
#'
#' # Example 3: non-stratified WRAT analysis with two time-to-event endpoints.
#' #   The "weights" method employed.
#' analysis.results <- WRAT(data = sample_data[,1:6],
#'                          method = "weights",
#'                          weights = c(1,1), calipers = 0.2,
#'                          clinical.threshold = c(0,0),
#'                          direction = c(1,1),
#'                          strata = NULL,
#'                          strata.weight = NULL)
#' cat("WR:",analysis.results$WR, "; P-value:", analysis.results$P,"\n")

WRAT <- function (data, method = c("weights","calipers"),
                  weights = NULL, calipers = NULL,
                  clinical.threshold = NULL,
                  direction = NULL,
                  strata = NULL, strata.weight = NULL) {
  ## check input information
  tryCatch({
    method <- match.arg(method)
  }, error = function(e) {
    stop("Argument 'method' must be 'weights' or 'calipers'", call. = FALSE)
  })

  if ((method == "calipers") & is.null(calipers)) {
    stop("Argument 'calipers' must be specified when using 'calipers' method")
  }

  ## strata and weights for strata
  if (is.null(strata)) strata = rep(1, nrow(data)) ## non-stratified analysis
  if (!is.null(strata.weight)) {
    oper.weight <- strata.weight
  } else {
    ## strata weights as number of participants if not specified
    stra <- unique(strata)
    oper.weight <- data.frame(table(strata))
    names(oper.weight) <- c('stra','weight')
  }
  oper.weight[order(oper.weight$stra),] # sort the stratum, increasing in index

  ## perform analysis
  num.outcome <- (ncol(data) - 2) / 2
  if (is.null(direction)) direction <- rep(1, num.outcome)
  # make larger number more favorable
  direction.matrix <- matrix(direction, nrow=nrow(data), ncol=num.outcome, byrow=T)
  data[,seq(3,(3+2*(num.outcome-1)),2)] <- data[,seq(3,(3+2*(num.outcome-1)),2)] *
    direction.matrix

  ## obtain pair.data.list
  pair.data.list <- list()
  for (s in 1:nrow(oper.weight)) {
    pair.data.list[[s]] <- data.to.pair(data[strata==oper.weight[s,1],])
  }
  ## combined pair.data for adaptive thresholds calculation
  pair.data.all <- do.call(rbind, pair.data.list)

  ## obtain adaptive thresholds
  ## weights method
  if (method == "weights") {
    if (is.null(weights)) {
      weights <- rep(1, num.outcome)
      message("Thresholds weights set to default: 1 for all outcomes")
    }

    if (is.null(calipers)) {
      calipers <- 0.2
      message("Thresholds caliper set to default: 20%")
    }

    adp.threshold <- adaptive.threshold.weights(pair.data = pair.data.all,
                                                num.outcome = num.outcome,
                                                weight = weights, caliper=calipers,
                                                clinical.threshold = clinical.threshold)

  ## calipers method
  } else {
    adp.threshold <- adaptive.threshold.calipers(pair.data = pair.data.all,
                                                 num.outcome = num.outcome,
                                                 calipers = calipers,
                                                 clinical.threshold = clinical.threshold)
  }

  ## perform pairwise comparisons for WRAT
  intermediate.res <- WRAT.comparison(data = data,
                                      pair.data.list = pair.data.list,
                                      threshold = adp.threshold$threshold,
                                      stage = adp.threshold$stage,
                                      direction=rep(1,num.outcome),
                                      strata = strata,
                                      strata.weight = oper.weight)


  res.output <- comparison.result(intermediate.res)

  return(res.output)
}
