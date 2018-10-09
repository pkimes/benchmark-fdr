#' Plot Mean Performance Metrics Across Replications
#'
#' Given a standardized metric table created using "standardize_results",
#' this function calculates the mean metric across replications and
#' generates a standard plot.
#' 
#' @param tsb standardized metric data.frame generated using
#'        standardize_results.
#' @param met name(s) of metric to plot - must be one of the performance
#'        metrics added by default with "addDefaultMetrics" or FWER or FPR.
#'        May be a vector of names of length 2 - in this case the first will
#'        be used for the x-axis and the second for the y-axis. 
#' @param filter_set character vector of "blabel" IDs of methods which
#'        shouldn't be included in the plot. Alternatively, a subset of
#'        methods can be chosen by filtering the data.table before passing
#'        to this function. (default = NULL)
#' @param merge_ihw logical whether IHW results should be merged across
#'        alpha values by matching "param.alpha" and "alpha" columns.
#'        (default = TRUE)
#' @param clean_names logical whether to clean-up method names in 'blabel'
#'        column. The tsb table must only contain one column with each of the
#'        following labels, else method labels will not be changed: 'ashq',
#'        'bh', 'bl', 'ihw', 'lfdr', 'qvalue', 'fdrreg-e',
#'        'fdrreg-t'. (default = FALSE)
#' @param errorBars logical indicating whether to include error bars 
#' @param palette data.frame containing the color palette - should contain 
#'        three columns: 'Method', 'col', and 'lty'
#' @param facetMethodType logical indicating whether to facet the plot into 
#'        two panels - methods that use only one piece of information, and 
#'        those that use more than just the p-value (default = FALSE)
#' @param diffplot logical indicating whether 'value' in tsb is a difference
#'        between informative and uninformative covariates. (default = FALSE)
#' @param grpVars vector of character names of additional columns of tsb to keep 
#'
#' @return
#' a ggplot object.
#' 
#' @author Patrick Kimes & Keegan Korthauer
plotsim_average <- function(tsb, met, filter_set = NULL, merge_ihw = TRUE,
                            clean_names = FALSE, errorBars=FALSE,
                            palette = candycols, facetMethodType = FALSE,
                            diffplot = FALSE, grpVars = NULL){
   if (length(met)>2)
     stop("Can only plot 2 metrics at a time")
  
    ## cacluate mean per replication
    if(!is.null(grpVars)){
      tsba <- tsb %>%
        group_by(blabel, performanceMetric, alpha, param.alpha, key, 
                 .dots = grpVars) %>%
        summarize(n = sum(!is.na(value)),
                  se = sd(value, na.rm = TRUE) / sqrt(n),
                  value = mean(value, na.rm = TRUE))
    }else{
      tsba <- tsb %>%
        group_by(blabel, performanceMetric, alpha, param.alpha, key) %>%
        summarize(n = sum(!is.na(value)),
                  se = sd(value, na.rm = TRUE) / sqrt(n),
                  value = mean(value, na.rm = TRUE))
   }
   
   tsba$value[tsba$n == 0] <- NA
   tsba$se[tsba$n == 0] <- NA
    
    ## remove methods if any specified
    if (!is.null(filter_set)) {
        tsba <- filter(tsba, !(blabel %in% filter_set))
    }

    ## group IHW methods if specified
    if (merge_ihw) {
        tsba <- filter(tsba, is.na(param.alpha) | (param.alpha == alpha))
        tsba$blabel[grepl("^ihw-", tsba$blabel)] <- "ihw"
    }
    
    # filter by performance metric
    tsba_m <- NULL 
    for (m in seq_along(met)){
      tmp <- tsba %>% 
        filter(performanceMetric == met[m]) %>%
        dplyr::mutate(Method = gsub("-df03", "", blabel)) 
      if (m > 1){
        tsba_m <- left_join(tsba_m, tmp, 
                            by = c("Method", "alpha", "n", 
                                   "param.alpha", "blabel"))
      }else{
        tsba_m <- tmp
      }
    }
    tsba <- tsba_m
           
    if (clean_names) {
        ulabs <- unique(tsba$blabel)
        vlabs <- c('ashq', 'bh', 'bl', 'ihw', 'lfdr', 'qvalue')
        clabs <- c("ASH q-value", "Benjamini-Hochberg", "Boca-Leek", "IHW", "local FDR", "Storey's q-value")
        if (any(grepl("fdrreg", ulabs))) {
            vlabs <- c(vlabs, 'fdrreg-e', 'fdrreg-t')
            clabs <- c(clabs, "FDRreg (emp)", "FDRreg (theor)")
        }
        if (any(grepl("adapt", ulabs))) {
            vlabs <- c(vlabs, 'adapt-glm', 'adapt-gam')
            clabs <- c(clabs, 'AdaPT (GLM)', 'AdaPT (GAM)')
        }
        lcnts <- sapply(vlabs, function(x) { grep(paste0("^", x), ulabs, value=TRUE) })
        if (!is(lcnts, "list")) {
            names(lcnts) <- clabs
            tsba$blabel <- do.call(forcats::fct_recode, c(list("f" = tsba$blabel), as.list(lcnts)))
        }
    }

        
    # add color palette
    tsba <- dplyr::left_join(tsba, palette, by="Method") 
  
    col <- as.character(tsba$col)
    names(col) <- as.character(tsba$Method)
    
    lty <- as.character(tsba$lty)
    names(lty) <- as.character(tsba$Method)
    
    # add method type
    tsba <- tsba %>%
      mutate(Type = ifelse(Method %in% c("unadjusted", "bonf", "bh", "qvalue"),
                           "Univariate (p-value only)", "Multivariate"))
    if(length(met)==1){
      gp <- tsba %>%
        ggplot(aes(x = alpha, y = value, color = Method)) +
        geom_line(alpha = 0.85, aes(linetype=Method)) +
        theme_classic() +
        theme(axis.title = element_text(face="bold"),
              plot.title = element_text(face="bold")) +
        expand_limits(x = 0) +
        scale_x_continuous("alpha cutoff", breaks=seq(0, 1, by=0.01)) +
        ylab(met) +
        ggtitle(ifelse(diffplot,
                       paste0("Mean Difference ", met, " Over ", max(tsba$n), " Replications"),
                       paste0("Mean ", met, " Over ", max(tsba$n), " Replications"))) +
        scale_color_manual(values = col) +
        scale_linetype_manual(values = lty)
    }else{
      gp <- tsba %>% 
        ggplot(aes(x = value.x, y = value.y, color = Method)) +
        geom_line(alpha = 0.85, aes(linetype=Method)) +
        theme_classic() +
        theme(axis.title = element_text(face="bold"),
              plot.title = element_text(face="bold")) +
        expand_limits(x = 0) +
        scale_x_continuous(breaks=seq(0, 1, by=0.01)) +
        ylab(met[2]) +
        xlab(met[1]) +
        ggtitle("Mean ROC Curve over 100 Replications") +
        scale_color_manual(values = col) +
        scale_linetype_manual(values = lty)
    }
    
    if(facetMethodType){
      gp <- gp +
        facet_wrap( ~ Type)
    }
        
    if(errorBars & length(met)==1){
      gp <- gp + geom_errorbar(width=0.0025, alpha=0.5,
                               aes(ymin=value-se, ymax=value+se))
    }

    ## use percentage on y-axis labels when appropriate
    if (sum(met %in% c("FDR", "FNR", "TPR", "TNR", "FWER", "rejectprop")) > 0){
      if (length(met)==1){
        gp <- gp +
            scale_y_continuous(ifelse(diffplot, paste(met, "(informative - uninformative)"), met),
                               labels=scales::percent)
      }else{
        gp <- gp +
          scale_x_continuous(met[1], labels=scales::percent) +
          scale_y_continuous(met[2], labels=scales::percent) 
      }
    } else if (diffplot) {
        gp <- gp + ylab(paste(met, "(informative - uninformative)")) 
    }

    
    if (diffplot) {
        gp <- gp + expand_limits(y = 0)
        gp <- gp + geom_hline(yintercept = 0, lty = 2, color = "blue", alpha = 1/2)
    } else {
        ## include 0% or 100% in plotting range depending on metric
        if (sum(met %in% c("FDR", "FNR", "FWER", "rejectprop")) > 0){
          if (length(met)==1)
            gp <- gp + expand_limits(y = 0)
        }
        if (sum(met %in% c("TNR")) > 0){
          if (length(met)==1)
            gp <- gp + expand_limits(y = 1)
        }

        ## add identity line for FPR/FDR plotting
        if (sum(met == "FDR")>0) {
          if (length(met)==1)
            gp <- gp +
                geom_abline(intercept = 0, slope = 1, lty = 2, color = "blue", alpha = 1/2)
        }
    }
    
    gp
}
