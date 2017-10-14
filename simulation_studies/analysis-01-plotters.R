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
#' a data.table of the following performance metrics:
#' 'rejections', 'TPR', 'TNR', 'FPR', 'FNR', 'FWER'.
#' 
#' @author Patrick Kimes
standardize_results <- function(res, alpha = seq(0.01, 0.10, 0.01)) {
    
    sbl <- lapply(res, addDefaultMetrics)
    sbl <- lapply(sbl, addPerformanceMetric,
                  evalMetric = "FWER", assay = "qvalue",
                  evalFunction =
                      function( query, truth, alpha = 0.1) {
                          any((query < alpha) & !as.logical(truth))
                      })
    sbl <- lapply(sbl, addPerformanceMetric,
                  evalMetric = "rejectprop", assay = "qvalue",
                  evalFunction =
                      function( query, truth, alpha = 0.1) {
                          mean(query < alpha)
                      })
    tsb <- lapply(sbl, estimatePerformanceMetrics,
                  alpha = alpha, tidy = TRUE)
    
    rbindlist(tsb, idcol = "rep")
}


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
#'
#' @return
#' a ggplot object.
#' 
#' @author Patrick Kimes
average_plots <- function(tsb, met, filter_set = NULL, merge_ihw = TRUE,
                          clean_names = FALSE) {

    ## cacluate mean over replications
    tsba <- tsb %>%
        group_by(blabel, performanceMetric, alpha, param.alpha, key) %>%
        summarize(value = mean(value), n = n())

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


#' Plot Adjusted P-Values Against Stratifying Variable
#'
#' Create histograms showing distribution of tests across a
#' stratifying variable (`strat`), colored by the adjusted 
#' p-values for each method. The plot is similar to Figure 1 from
#' the ASH paper, used to show the relationship between the adjusted
#' p-values and the effect size of a contrast being tested.
#' 
#' @param sba data.frame with qvalue assay from a
#'        single SummarizedBenchmark object.
#' @param df data.frame of simulation data used to generate `sb`.
#' @param strat character string corresponding to column in `df` to
#'        be used as stratifying variable over which adjusted p-values
#'        should be plotted. (default = "effect_size")
#' @param clean_names logical whether to clean-up method names. The sba table
#'        must only contain one column for each of the following labels,
#'        else method labels will not be changed: 'ashs',
#'        'bh', 'bl', 'ihw', 'lfdr', 'qvalue', 'scott-empirical',
#'        'scott-theoretical'. (default = FALSE)
#'
#' @return
#' a ggplot object
#' 
#' @author Patrick Kimes
single_plots <- function(sba, df, strat = "effect_size", clean_names = FALSE) {

    ## check compatibility of tables
    stopifnot(is.data.frame(sba))
    stopifnot(nrow(sba) == nrow(df))

    ## verify stratifying variable is valid column
    stopifnot(strat %in% names(df))

    ## merge adjusted p-values with stratifying variable
    sba <- cbind(sba, strat = df[, strat])
    names(sba)[names(sba) == "strat"] <- strat
    sba <- melt(sba, id.vars = strat, variable.name = "method",
                value.name = "adjp")

    ## relabel methods with nicer names if possible
    if (clean_names) {
        ulabs <- unique(sba$method)
        vlabs <- c('ashs', 'bh', 'bl', 'ihw', 'lfdr',
                   'qvalue', 'scott.empirical', 'scott.theoretical')
        clabs <- c("ASH s-value", "Benjamini-Hochberg", "Boca-Leek", "IHW", "local FDR",
                   "Storey's q-value", "FDRreg (emp)", "FDRreg (theor)")
        lcnts <- sapply(vlabs, function(x) { grep(paste0("^", x), ulabs, value=TRUE) })
        if (!is(lcnts, "list")) {
            names(lcnts) <- clabs
            sba$method <- do.call(forcats::fct_recode, c(list("f" = sba$method), as.list(lcnts)))
        }
    }

    ## check for negative adjusted p-values and throw warning; still plot with 0s
    if (any(sba$adjp < 0)) {
        negp <- which(sba$adjp < 0)
        nneg <- length(negp)
        mneg <- paste(unique(sba$method[negp]), collapse=",")
        warning("adjusted p-values include ", nneg, " negative values ",
                "from methods: ", mneg, "! ",
                "setting these to 0 for plotting.")
        sba$adjp[sba$adjp < 0] <- 0
    }

    ## plot adjusted p-values against stratifying variable in df
    sba %>%
        ggplot(aes(x = get(strat), fill = .p2disc(adjp))) +
        geom_histogram(color='black', size=1/4, bins=40,
                       boundary=0, position="stack") +
        xlab(strat) + 
        viridis::scale_fill_viridis("adjusted\np-value", discrete=TRUE,
                                    option="plasma", direction=-1,
                                    begin=.1, end=.9, drop=FALSE) +
        theme_classic() +
        theme(axis.title = element_text(face="bold"),
              plot.title = element_text(face="bold")) +
        ggtitle("Adjusted P-Values") +
        facet_wrap(~ method)
}


## helper; bin adjusted p-vals, flip order of bins
.p2disc <- function(x) {
    cut(x, breaks=c(0, .01, .05, .10, .20, 1), include.lowest=TRUE) %>%
        factor(., levels=rev(levels(.)))
}


