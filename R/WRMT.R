#' Organize dataset for pairwise comparison
#'
#' @param data Input data frame
#'
#' @returns A data frame for pairwise comparison
#' @export
#'
#' @examples
data.to.pair <- function (data) {
  left.part <- dplyr::slice(data, rep(1:dplyr::n(), each = nrow(data)))
  right.part <- dplyr::slice(data, rep(1:dplyr::n(), times = nrow(data)))
  output <- cbind(left.part,right.part)
  return(output)
}


#' Comparison for a single stage with potential right censoring
#'
#' @param t1 A numeric vector. Observed outcome for target participants.
#' @param t2 A numeric vector. Observed outcome for participants for comparison.
#' @param c1 A numeric vector. Censoring indicator for target participants.
#' @param c2 A numeric vector. Censoring indicator for participants for comparison.
#' @param threshold A numeric value. Threshold used for comparison at this stage,
#' differences smaller than this threshold are considered undetermined.
#'
#' @returns A numeric vector that indicates the comparison results for this stage:
#' 1 for wins by target participants, -1 for losses by target participants,
#' and 0 for undetermined pairs.
#' @export
#'
#' @examples
WRMT.one.stage.sur <- function (t1,t2,c1,c2,threshold) {
  c00 <- (c1+c2) == 0 ## both not censored
  threshold.indi <- abs(t1-t2) >= threshold
  direction.indi <- sign(t1 - t2)

  c10 <- c1 > c2 ## c1 censored, c2 not censored
  c10.dir <- (t1-t2)>=threshold

  c01 <- c1 < c2 ## c2 censored, c1 not censored
  c01.dir <- (t2-t1)>=threshold

  output <- c00 * threshold.indi * direction.indi + c10 * c10.dir - c01 * c01.dir
  return(output)
}

#' Pairwise comparisons for WRAT
#'
#' @description
#' This function performs pairwise comparisons for WRMT analysis. To be called
#' by function \code{WRMT()}.
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
#'
#' @returns A list of comparison results:
#' \describe{
#'   \item{pair.data.list}{A list that contains data frames of pairwise comparison
#'   results and data used for comparison. Each stratum corresponds to one data frame.}
#'   \item{stage.results.list}{A list that contains matrices of comparison results of each
#'   stage. Each stratum corresponds to one matrix.
#'   1 for win, -1 for loss, 0 for tie, and NA for results determined by previous stages.}
#'   \item{treatment}{A numeric vector indicates the treatment status.}
#'   \item{strata}{A numeric vector indicates the stratum that each participant belongs to.}
#'   \item{strata.weight}{A data frame contains weighting information
#'   for stratification analysis.}
#' }
#' @export
#'
#' @examples
#' # Example 1: stratified comparison
#' strata.compare <- WRMT.comparison(data = sample_data[,1:6],
#'                                   threshold = c(0,0), stage = c(1,2),
#'                                   direction = c(1,1),
#'                                   strata = sample_data$strata,
#'                                   strata.weight = NULL)
#'
#' # Example 2: non-stratified comparison
#' nonstrata.compare <- WRMT.comparison(data = sample_data[,1:6],
#'                                      threshold = c(0,0), stage = c(1,2),
#'                                      direction = c(1,1),
#'                                      strata.weight = NULL)
#'
#'
WRMT.comparison <- function (data, threshold, stage=NULL, direction=NULL,
                             strata = NULL, strata.weight = NULL) {

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

  ## get endpoint information (number of endpoints, favorable values, etc.)
  num.outcome <- (ncol(data) - 2) / 2
  if (is.null(direction)) direction <- rep(1, num.outcome)
  # make larger number more favorable
  direction.matrix <- matrix(direction, nrow=nrow(data), ncol=num.outcome, byrow=T)
  data[,seq(3,(3+2*(num.outcome-1)),2)] <- data[,seq(3,(3+2*(num.outcome-1)),2)] *
    direction.matrix

  ## check stage information correctly specified
  if (is.null(stage) & (length(threshold)%%num.outcome == 0)) {
    stage <- rep(seq(1,num.outcome),length(threshold)/num.outcome)
  } else if (is.null(stage)) {
    stop('Stage needs to be specified under current setting')
  }

  ## pairwise comparison for each stratum
  pair.data.list <- list()
  stage.results.list <- list()
  for (s in 1:nrow(oper.weight)) {
    current.strata.data <- data[strata==oper.weight[s,1],]
    pair.data <- data.to.pair(current.strata.data)
    pair.data$result <- 0 # comparison results
    stage.results <- matrix(data=NA, nrow = nrow(pair.data), ncol = length(stage))

    col.num <- 2 + 2*num.outcome

    for (i in 1:length(threshold)) {
      current.sub <- pair.data$result == 0
      t1.index <- 1 + 2 * stage[i]
      t2.index <- col.num + t1.index
      c1.index <- 2 + 2 * stage[i];
      c2.index <- col.num + c1.index
      current.stage <- WRMT.one.stage.sur(t1=pair.data[current.sub,t1.index],
                                          t2=pair.data[current.sub,t2.index],
                                          c1=pair.data[current.sub,c1.index],
                                          c2=pair.data[current.sub,c2.index],
                                          threshold=threshold[i])
      pair.data$result[current.sub] <- current.stage
      stage.results[current.sub,i] <- current.stage
    }

    pair.data.list[[s]] <- pair.data
    stage.results.list[[s]] <- stage.results
    names(pair.data.list)[s] <- paste0("stratum", oper.weight[s,1])
    names(stage.results.list)[s] <- paste0("stratum", oper.weight[s,1])

  }

  ## organize output list
  output <- list(pair.data.list = pair.data.list,
                 stage.results.list = stage.results.list,
                 treatment = data[,2],
                 strata = strata,
                 strata.weight = oper.weight)

  return(output)
}


#' Summarize comparison results for outputs
#'
#' @description
#' This function analyzes pairwise comparisons results provided by
#' \code{WRMT.comparison()} or \code{WRAT.comparison()}. To be called
#' by function \code{WRMT()} or \code{WRAT()}.
#'
#' @param comparison.output
#' A list returned by \code{WRMT.comparison()} or \code{WRAT.comparison()}.
#' Must contain elements \code{pair.data}, \code{stage.results},
#' \code{treatment}, and \code{strata}.
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
#'
#' @export
#'
#' @examples
#' # Example 1: stratified comparison
#' strata.compare <- WRMT.comparison(data = sample_data[,1:6],
#'                                   threshold = c(0,0), stage = c(1,2),
#'                                   direction = c(1,1),
#'                                   strata = sample_data$strata,
#'                                   strata.weight = NULL)
#' analysis.result <- comparison.result(strata.compare)
#' cat("WR:",analysis.result$WR, "; P-value:", analysis.result$P,"\n")
#'
#' # Example 2: non-stratified comparison
#' nonstrata.compare <- WRMT.comparison(data = sample_data[,1:6],
#'                                      threshold = c(0,0), stage = c(1,2),
#'                                      direction = c(1,1),
#'                                      strata = NULL,
#'                                      strata.weight = NULL)
#' analysis.result <- comparison.result(nonstrata.compare)
#' cat("WR:",analysis.result$WR, "; P-value:", analysis.result$P,"\n")
#'
comparison.result <- function (comparison.output) {

  ## extract strata information
  strata <- comparison.output$strata
  oper.weight <- comparison.output$strata.weight

  ## stratified analysis for comparison results within each stratum
  strata.WR <- strata.stat <- strata.var <- rep(0,nrow(oper.weight))
  strata.win <- strata.loss <- rep(0,nrow(oper.weight))
  strata.stage.output <- list() # list for stage level decomposition

  for (s in 1:nrow(oper.weight)) {
    # extract pair.data for current stratum
    pair.data.input <- comparison.output[["pair.data.list"]][[s]]
    org.treat <- with(comparison.output,
                      treatment[strata == oper.weight[s,1]])

    col.num <- (ncol(pair.data.input)-1)/2
    win <- sum((pair.data.input[,2]==1) * (pair.data.input[,(2+col.num)]==0) *
                 (pair.data.input[,(1 + 2*col.num)]==1))
    loss <- sum((pair.data.input[,2]==1) * (pair.data.input[,(2+col.num)]==0) *
                  (pair.data.input[,(1 + 2*col.num)]== -1))
    treat.sum <- sum((pair.data.input[,2]==1) * pair.data.input[,(1 + 2*col.num)])
    WR <- win/loss
    var.data <- pair.data.input[,c(1,(1 + 2*col.num))]
    var.data.oper <- aggregate(var.data,by=list(var.data[,1]),FUN=sum)[,3]
    var <- sum(var.data.oper^2)*sum(org.treat==1)*sum(org.treat==0)/
      (length(org.treat)*(length(org.treat)-1))
    strata.win[s] <- win
    strata.loss[s] <- loss
    strata.WR[s] <- WR
    strata.stat[s] <- treat.sum
    strata.var[s] <- var

    ## stage level wins, losses, and ties
    # extract stage.results for stage-level information
    stage.results.input <- comparison.output[["stage.results.list"]][[s]]
    num.stage <- ncol(stage.results.input)
    stage.results.input <- cbind(pair.data.input[,c(2,(2+col.num))],stage.results.input)
    stage.output <- matrix(NA, nrow=num.stage, ncol=3)
    num.between.group <- sum((stage.results.input[,1]==1) * (stage.results.input[,2]==0))
    for (stage.index in 1:num.stage) {
      win <- sum((stage.results.input[,1]==1) * (stage.results.input[,2]==0) *
                   (stage.results.input[,(2 + stage.index)]==1), na.rm = T)
      loss <- sum((stage.results.input[,1]==1) * (stage.results.input[,2]==0) *
                    (stage.results.input[,(2 + stage.index)]==-1), na.rm = T)
      tie <- sum((stage.results.input[,1]==1) * (stage.results.input[,2]==0) *
                   (stage.results.input[,(2 + stage.index)]==0), na.rm = T)
      stage.output[stage.index,] <- c(win, tie, loss) / num.between.group
    }
    stage.output <- as.data.frame(stage.output) * 100
    names(stage.output) <- c("wins(%)","ties(%)","losses(%)")
    strata.stage.output[[s]] <- stage.output
    names(strata.stage.output)[s] <- paste0("stratum",oper.weight[s,1])
  }

  ## integrate information from strata
  WR <- sum((strata.win/oper.weight[,2])) / sum((strata.loss/oper.weight[,2]))
  # WR_old <- weighted.mean(strata.WR,w=wr.weight)
  stat <- sum(strata.stat)
  var <- sum(strata.var)
  P <- 2*pnorm(abs(stat/sqrt(var)),lower.tail = F)

  strata.summary <- data.frame(strata = oper.weight[,1],
                               weight = oper.weight[,2]/sum(oper.weight[,2]),
                               unadj.WR = strata.WR
                               )

  ## simplify output if only one stratum
  if (nrow(oper.weight) == 1) {
    strata.stage.output <- strata.stage.output[[1]]
    strata.summary <- NULL
  }

  output <- list(WR=WR, P=P,
                 Stage.decompose = strata.stage.output,
                 strata.summary = strata.summary,
                 FS.stat=stat, FS.stat.var=var)
  return(output)
}


#' Obtain adaptive thresholds with universal caliper c and weights w
#' #' @description
#' This function calculates adaptive thresholds for WRAT analysis with universal
#' caliper c and weights w. To be called by function \code{WRAT()}.
#'
#' @param pair.data A data frame for pairwise comparison,
#' returned by \code{data.to.pair()}.
#' @param num.outcome A numeric value of the number of outcomes.
#' @param weight A numeric vector indicates the weights w for adaptive thresholds.
#' @param caliper A numeric value indicates the universal caliper c for adaptive thresholds.
#' @param clinical.threshold
#' A numeric vector represents the minimal clinical thresholds.
#' Default to 0 for all endpoints.
#'
#' @returns A list of outputs:
#' \describe{
#'   \item{threshold}{A numeric vector of obtained adaptive thresholds.}
#'   \item{stage}{A numeric vector indicating outcome for each comparison stage.}
#' }
#' @export
#'
#' @examples
adaptive.threshold.weights <- function (pair.data, num.outcome,
                                        weight, caliper=0.2,
                                        clinical.threshold = NULL) {
  col.num <- ncol(pair.data) / 2
  adp.thresholds <- NULL
  ## calculate adaptive thresholds for endpoints
  for (j in 1:num.outcome) {
    index <- 3+2*(j-1)
    diff <- pair.data[,index] - pair.data[,(index+col.num)]
    diff.df <- diff[diff > 0] ## non-zero difference only
    diff.quantile <- quantile(diff.df, probs = caliper)
    adp.thresholds[j] <- diff.quantile/weight[j]
  }

  ## involve clinical thresholds
  if (is.null(clinical.threshold)) clinical.threshold <- rep(0,num.outcome)
  adp.thresholds <- pmax(adp.thresholds, clinical.threshold)
  adp.thresholds[seq(num.outcome+1,2*num.outcome)] <- clinical.threshold

  stage <- rep(seq(1,num.outcome),2)

  output <- list(threshold = adp.thresholds,
                 stage = stage)
  return(output)
}


#' Obtain adaptive thresholds with caliper c specified for each outcome
#'
#' @description
#' This function calculates adaptive thresholds for WRAT analysis with
#' caliper c specified for each outcome. To be called by function \code{WRAT()}.
#'
#' @param pair.data A data frame for pairwise comparison,
#' returned by \code{data.to.pair()}.
#' @param num.outcome A numeric value of the number of outcomes.
#' @param caliper A numeric vector indicates the caliper c used for each outcome.
#' @param clinical.threshold
#' A numeric vector represents the minimal clinical thresholds.
#' Default to 0 for all endpoints.
#'
#' @returns A list of outputs:
#' \describe{
#'   \item{threshold}{A numeric vector of obtained adaptive thresholds.}
#'   \item{stage}{A numeric vector indicating outcome for each comparison stage.}
#' }
#' @export
#'
#' @examples
adaptive.threshold.calipers <- function (pair.data, num.outcome, calipers,
                                         clinical.threshold = NULL) {
  col.num <- ncol(pair.data) / 2
  adp.thresholds <- NULL
  ## calculate adaptive thresholds for endpoints
  for (j in 1:num.outcome) {
    index <- 3+2*(j-1)
    index.c <- index + 1
    diff <- pair.data[,index] - pair.data[,(index+col.num)]
    diff.df <- diff[diff > 0] ## non-zero difference only
    diff.quantile <- quantile(diff.df, probs = calipers[j])
    adp.thresholds[j] <- diff.quantile
  }

  ## involve clinical thresholds
  if (is.null(clinical.threshold)) clinical.threshold <- rep(0,num.outcome)
  adp.thresholds <- pmax(adp.thresholds, clinical.threshold)
  adp.thresholds[seq(num.outcome+1,2*num.outcome)] <- clinical.threshold

  stage <- rep(seq(1,num.outcome),2)

  output <- list(threshold = adp.thresholds,
                 stage = stage)
  return(output)
}


#' Pairwise comparisons for WRAT
#'
#' @description
#' This function performs pairwise comparisons for WRAT analysis. To be called
#' by function \code{WRAT()}.
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
#' @param pair.data.list A list of data frames organized for
#' stratified pairwise comparisons.
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
#'
#' @returns A list of comparison results:
#' \describe{
#'   \item{pair.data.list}{A list that contains data frames of pairwise comparison
#'   results and data used for comparison. Each stratum corresponds to one data frame.}
#'   \item{stage.results.list}{A list that contains matrices of comparison results of each
#'   stage. Each stratum corresponds to one matrix.
#'   1 for win, -1 for loss, 0 for tie, and NA for results determined by previous stages.}
#'   \item{treatment}{A numeric vector indicates the treatment status.}
#'   \item{strata}{A numeric vector indicates the stratum that each participant belongs to.}
#'   \item{strata.weight}{A data frame contains weighting information
#'   for stratification analysis.}
#' }
#' @export
#'
#' @examples
#'
#'
WRAT.comparison <- function (data, pair.data.list, threshold, stage=NULL, direction=NULL,
                             strata = NULL, strata.weight = NULL) {

  ## strata and weights for strata
  oper.weight <- strata.weight

  ## get endpoint information (number of endpoints, favorable values, etc.)
  num.outcome <- (ncol(data) - 2) / 2

  ## pairwise comparison for each stratum
  stage.results.list <- list()
  for (s in 1:nrow(oper.weight)) {
    pair.data <- pair.data.list[[s]]
    pair.data$result <- 0 # comparison results
    stage.results <- matrix(data=NA, nrow = nrow(pair.data), ncol = length(stage))

    col.num <- 2 + 2*num.outcome

    for (i in 1:length(threshold)) {
      current.sub <- pair.data$result == 0
      t1.index <- 1 + 2 * stage[i]
      t2.index <- col.num + t1.index
      c1.index <- 2 + 2 * stage[i];
      c2.index <- col.num + c1.index
      current.stage <- WRMT.one.stage.sur(t1=pair.data[current.sub,t1.index],
                                          t2=pair.data[current.sub,t2.index],
                                          c1=pair.data[current.sub,c1.index],
                                          c2=pair.data[current.sub,c2.index],
                                          threshold=threshold[i])
      pair.data$result[current.sub] <- current.stage
      stage.results[current.sub,i] <- current.stage
    }

    pair.data.list[[s]] <- pair.data
    stage.results.list[[s]] <- stage.results
    names(pair.data.list)[s] <- paste0("stratum", oper.weight[s,1])
    names(stage.results.list)[s] <- paste0("stratum", oper.weight[s,1])

  }

  ## organize output list
  output <- list(pair.data.list = pair.data.list,
                 stage.results.list = stage.results.list,
                 treatment = data[,2],
                 strata = strata,
                 strata.weight = oper.weight)

  return(output)
}
