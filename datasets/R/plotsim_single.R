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
#'        else method labels will not be changed: 'ashq',
#'        'bh', 'bl', 'ihw', 'lfdr', 'qvalue', 'fdrreg-e',
#'        'fdrreg-t'. (default = FALSE)
#'
#' @return
#' a ggplot object.
#' 
#' @author Patrick Kimes
plotsim_single <- function(sba, df, strat = "effect_size", clean_names = FALSE) {

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
        vlabs <- c('ashq', 'bh', 'bl', 'ihw', 'lfdr', 'qvalue')
        clabs <- c("ASH q-value", "Benjamini-Hochberg", "Boca-Leek", "IHW", "local FDR", "Storey's q-value")
        if (any(grepl("fdrreg", ulabs))) {
            vlabs <- c(vlabs, 'fdrreg-e', 'fdrreg-t')
            clabs <- c(clabs, "FDRreg (emp)", "FDRreg (theor)")
        }
        lcnts <- sapply(vlabs, function(x) { grep(paste0("^", x), ulabs, value=TRUE) })
        if (!is(lcnts, "list")) {
            names(lcnts) <- clabs
            sba$method <- do.call(forcats::fct_recode, c(list("f" = sba$method), as.list(lcnts)))
        }
    }

    ## check for negative adjusted p-values and throw warning; still plot with 0s
    if (any(sba$adjp < 0, na.rm=TRUE)) {
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


