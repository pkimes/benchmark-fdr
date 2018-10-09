#' Standardize Simulation Results
#'
#' Given a list of SummarizedBenchmark objects with a "qvalue" assay,
#' this function calculates performance metrics and returns a single
#' data.table that can be used for plotting with the average_plots function.
#' 
#' @param res list of SummarizedBenchmark objects with "qvalue" assay that
#'        correspond to repilications of the same simulation setting. 
#' @param alpha sequence of alpha cutoffs at which metrics should be
#'        calculated. (default = seq(0.01, 0.10, 0.01))
#' 
#' @return
#' a tibble of the following performance metrics:
#' 'rejections', 'TPR', 'TNR', 'FPR', 'FNR', 'FWER'.
#' 
#' @author Patrick Kimes & Keegan Korthauer
plotsim_standardize <- function(res, alpha = seq(0.01, 0.10, 0.01)) {
    
    sbl <- lapply(res, addDefaultMetrics)
    sbl <- lapply(sbl, addPerformanceMetric,
                  evalMetric = "FWER", assay = "qvalue",
                  evalFunction =
                      function( query, truth, alpha = 0.1) {
                          any((query < alpha) & !as.logical(truth), na.rm = TRUE)
                      })
    sbl <- lapply(sbl, addPerformanceMetric,
                  evalMetric = "rejectprop", assay = "qvalue",
                  evalFunction =
                      function( query, truth, alpha = 0.1) {
                          mean(query < alpha, na.rm = TRUE)
                      })
    
    # add FPR metric for ROC curve
    sbl <- lapply(sbl, addPerformanceMetric,
                  assay="qvalue", evalMetric="FPR",
                  evalFunction = 
                    function( query, truth, alpha=0.1 ){
                    sum( query <= alpha & truth == 0, na.rm = TRUE ) / 
                         sum( truth == 0, na.rm = TRUE )
                    })
    
    # override default FDR perf. metric st fdr is zero (instead of NA) when
    # there are zero rejections
    sbl <- lapply(sbl, addPerformanceMetric,
                  evalMetric = "FDR", assay = "qvalue",
                  evalFunction =
                    function( query, truth, alpha = 0.1) {
                     rejections <- sum( query <= alpha, na.rm = TRUE )
                     fdr <- sum( query <= alpha & truth == 0, na.rm = TRUE ) / 
                        sum( query <= alpha, na.rm = TRUE )
                     fdr[rejections == 0] <- 0
                     return(fdr)
                    })
    
    tsb <- lapply(sbl, estimatePerformanceMetrics,
                  alpha = alpha, tidy = TRUE)
    
    tsb <- bind_rows(tsb, .id = "rep")
    as_tibble(tsb)
}
