#' Plot Mean Performance Metrics Across Replications
#'
#' Given a standardized metric table created using "standardize_results",
#' this function calculates the mean metric across replications and
#' generates a standard plot.
#' 
#' @param tsb standardized metric data.table generated using
#'        standardize_results.
#' @param met name of metric to plot - must be one of the performance
#'        metrics added by default with "addDefaultMetrics" or FWER.
#' @param filter_set character vector of "blabel" IDs of methods which
#'        shouldn't be included in the plot. Alternatively, a subset of
#'        methods can be chosen by filtering the data.table before passing
#'        to this function. (default = NULL)
#' @param merge_ihw logical whether IHW results should be merged across
#'        alpha values by matching "param.alpha" and "alpha" columns.
#'        (default = TRUE)
#' @param clean_names logical whether to clean-up method names in 'blabel'
#'        column. The tsb table must only contain one column with each of the
#'        following labels, else method labels will not be changed: 'ashs',
#'        'bh', 'bl', 'ihw', 'lfdr', 'qvalue', 'scott-empirical',
#'        'scott-theoretical'. (default = FALSE)
#' @param errorBars logical indicating whether to include error bars 
#'
#' @return
#' a ggplot object.
#' 
#' @author Patrick Kimes
plotsim_average <- function(tsb, met, filter_set = NULL, merge_ihw = TRUE,
                            clean_names = FALSE, errorBars=FALSE) {

    ## cacluate mean over replications
    tsba <- tsb %>%
        group_by(blabel, performanceMetric, alpha, param.alpha, key) %>%
        summarize(n = n(),
                  se = sd(value)/sqrt(n),
                  value = mean(value))

    ## remove methods if any specified
    if (!is.null(filter_set)) {
        tsba <- filter(tsba, !(blabel %in% filter_set))
    }

    ## group IHW methods if specified
    if (merge_ihw) {
        tsba <- filter(tsba, is.na(param.alpha) | (param.alpha == alpha))
        tsba$blabel[grepl("^ihw-", tsba$blabel)] <- "ihw"
    }

    if (clean_names) {
        ulabs <- unique(tsba$blabel)
        vlabs <- c('ashs', 'bh', 'bl', 'ihw', 'lfdr',
                   'qvalue', 'scott-empirical', 'scott-theoretical')
        clabs <- c("ASH s-value", "Benjamini-Hochberg", "Boca-Leek", "IHW", "local FDR",
                   "Storey's q-value", "FDRreg (emp)", "FDRreg (theor)")
        lcnts <- sapply(vlabs, function(x) { grep(paste0("^", x), ulabs, value=TRUE) })
        if (!is(lcnts, "list")) {
            names(lcnts) <- clabs
            tsba$blabel <- do.call(forcats::fct_recode, c(list("f" = tsba$blabel), as.list(lcnts)))
        }
    }
    
    gp <- tsba %>%
        filter(performanceMetric == met) %>%
        ggplot(aes(x = alpha, y = value, color = blabel)) +
        geom_line(alpha = 1/3) +
        geom_point(alpha = 1) +
        theme_classic() +
        theme(axis.title = element_text(face="bold"),
              plot.title = element_text(face="bold")) +
        expand_limits(x = 0) +
        viridis::scale_color_viridis("Method", discrete = TRUE,
                                     guide = guide_legend(ncol = 2)) +
        scale_x_continuous("alpha cutoff", breaks=seq(0, 1, by=0.01)) +
        ylab(met) +
        ggtitle(paste0("Mean ", met, " Over ", max(tsba$n), " Replications")) +
        theme(legend.position = c(0.02, 0.98),
              legend.justification = c("left", "top"),
              legend.background = element_rect(fill=scales::alpha("gray90", 1/3),
                                               color="black"))
    if(errorBars){
      gp <- gp + geom_errorbar(width=0.0025, alpha=1/3,
                               aes(ymin=value-se, ymax=value+se))
    }

    ## use percentage on y-axis labels when appropriate
    if (met %in% c("FPR", "FNR", "TPR", "TNR", "FWER", "rejectprop")) {
        gp <- gp +
            scale_y_continuous(met, labels=scales::percent)
    }

    ## include 0% or 100% in plotting range depending on metric
    if (met %in% c("FPR", "FNR", "FWER", "rejectprop")) {
        gp <- gp + expand_limits(y = 0)
    }
    if (met %in% c("TPR", "TNR")) {
        gp <- gp + expand_limits(y = 1)
    }

    ## add identity line for FPR/FDR plotting
    if (met == "FPR") {
        gp <- gp +
            geom_abline(intercept = 0, slope = 1, lty = 2, color = "blue", alpha = 1/2)
    }

    gp
}
