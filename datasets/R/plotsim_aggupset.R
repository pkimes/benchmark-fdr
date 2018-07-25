#' Aggregated UpSet Plot
#'
#' Generate an UpSet plot showing the (rounded) median overlap between methods
#' across a collection of SummarizedBenchmark objects, e.g. corresponding
#' to simulation replicates.
#' 
#' @param res list of SummarizedBenchmark objects to be combined in the plot.
#' @param alpha significance threshold to use for distinguishing significant
#'        and non-significant tests.
#' @param supplementary logical whether plot is for supplementary materials.
#'        (default = FALSE)
#' @param return_list logical whether frequency list should be returned instead
#'        of upset plot. The returned list can be used to generate the upset plot
#'        using `upset(fromExpression(freq_list)`. This can be useful if the user
#'        wants to experiment with upset plot styles. (default = FALSE)
#' @param nintersects scalar value representing number of sets to look at.  
#'        Default is 40 (same as default in UpSetR package).
#' 
#' @return
#' an upset plot if `return_list = FALSE` (default), else a list of frequencies that
#' can be used to generate the same upset plot.
#'
#' @details
#' Note: this can get incredibly slow if the number of methods being compared is
#' large since the possible number of overlaps grows exponentially. Anecdotally,
#' in simulations comparing 9 methods (+ truth) with 20,000 tests takes approximately
#' 15 to 20 seconds for 20 replications, and 40 to 45 seconds with 100 replications.  
#' 
#' @import dplyr magrittr
#' @author Patrick Kimes
aggupset <- function(res, alpha, supplementary = FALSE, return_list = FALSE,
                     nintersects = 40) { 

    ## find significant hits at alpha cutoff for all replicates
    hits_tabs <- lapply(res, sb2hits, a = alpha, s = supplementary)

    ## replace NAs with 0s (not called significant)
    fails <- lapply(hits_tabs, sapply, function(x) { all(is.na(x)) })
    hits_tabs <- lapply(hits_tabs, function(x) { x[is.na(x)] <- 0; x })
    
    ## count up frequencies in each intersection
    n_cols <- unique(sapply(hits_tabs, ncol))
    if (length(n_cols) > 1) {
        stop("not all SummarizedBenchmarks have the same set of methods")
    }
    freq_tabs <- lapply(hits_tabs, hits2freq, nm = n_cols)

    ## convert anything that failed completely to NAs
    freq_tabs <- mapply(function(x, y) {
        if (!any(x)) { return(y) }
        failid <- make.names(names(x))[x]
        failid <- match(failid, names(y))
        y$freq[rowSums(y[, failid]) > 0] <- NA
        y
    }, x = fails, y = freq_tabs, SIMPLIFY = FALSE) 
    
    ## merge all freqs into single table - first rename 'freq' columns to 'freq.i' (i = 1..100)
    method_names <- setdiff(names(freq_tabs[[1]]), "freq")
    freq_tab <- mapply(function(itab, idx) { dplyr::rename(itab, !!(paste0("freq.", idx)) := freq) },
                       itab = freq_tabs, idx = 1:length(freq_tabs),
                       SIMPLIFY = FALSE) %>%
        purrr::reduce(dplyr::left_join, by = method_names)

    ## summarize across 100 replications of each setting
    freq_tab <- freq_tab %>%
        gather(repl, cnt, starts_with("freq")) %>%
        group_by_at(method_names) %>%
        summarize(freq_mean = round(mean(cnt, na.rm = TRUE))) %>%
        ungroup() %>%
        mutate(freq_mean = ifelse(is.nan(freq_mean), 0, freq_mean)) 
    
    ## convert binary design matrix to UpSetR format (method names separated by "&")
    freq_tab <- freq_tab %>%
        unite("design", method_names, sep = "&", remove = FALSE) %>%
        gather(method, val, -design, -freq_mean) %>%
        mutate(val = ifelse(val, method, "")) %>%
        spread(method, val) %>%
        select(-design) %>%
        unite("setname", method_names, sep = "&") %>%
        mutate(setname = setname %>% gsub("&{2,}", "&", .) %>%
                   gsub("^&", "", .) %>%
                   gsub("&$", "", .)) 

    ## convert to vector to pass to UpSetR package
    freq_list <- freq_tab$freq_mean
    names(freq_list) <- freq_tab$setname

    ## return frequency list if requested
    if (return_list) {
        return(freq_list)
    }
    
    ## draw upset plot if frequency list not returned
    upset(fromExpression(freq_list),
          nsets = n_cols,
          nintersects = nintersects,
          mb.ratio = c(0.55, 0.45),
          order.by = "freq",
          decreasing = TRUE,
          set.metadata = list(data = data.frame(sets = method_names,
                                                isTruth = grepl("truth", method_names)),
                              plots = list(
                                  list(type = "matrix_rows", 
                                       column = "isTruth",
                                       colors = c("TRUE" = "blue", "FALSE" = "gray"), 
                                       alpha = 0.2))))
}


#' Helper to Parse Significant Hits for Specified Alpha
#'
#' Determines which tests are significant for each method based
#' on a specified alpha threshold and returns as a binary data.frame
#' for easier downstream parsing.
#' 
#' @param x SummarizedBenchmark w/ qvalue assay.
#' @param a alpha cutoff.
#' @param s logical whether for supplementary materials or not.
#'
#' @return
#' data.frame of 0/1s; rows are test, columns are methods.
#'
#' @import dplyr magrittr
#' @author Patrick Kimes
sb2hits <- function(x, a, s) {
    ## make quick table of significant tests w/ groundTruth
    ht <- as_tibble(cbind((assay(x, "qvalue") < a) + 0,
                          truth = rowData(x)$qvalue))
    ## keep only IHW matching "alpha" parameter
    ihw_keep <- paste0("ihw-a", sprintf("%02i", 100 * a ))
    if (ihw_keep %in% names(ht)) {
        ## note - using mutate instead of rename so next 'select' call to drop
        ## extra "ihw-*" columns doesn't throw an error if correct alpha was only alpha
        ht <- dplyr::mutate_(ht, ihw = paste0("`", ihw_keep, "`"))
    }
    ht <- dplyr::select(ht, -dplyr::contains("ihw-"))
    ## if not plotting for supplementary materials, remove BL w/ multiple DoF 
    if (!s) {
        suppressWarnings({
            ht <- ht %>%
                dplyr::select(-one_of("bl-df02", "bl-df04", "bl-df05")) %>%
                dplyr::rename(bl = `bl-df03`)
        })
    }
    suppressWarnings({
        ht <- dplyr::select(ht, -one_of("unadjusted"))
    })
    as.data.frame(ht)
}


#' Helper to Count Intersection Frequencies Across Methods
#'
#' Counts overlaps/intersections between methods based on the
#' binary data.frame generated using the `sb2hits()` function.
#' This function is just a wrapper to the `Counter()` function
#' in the `UpSetR` package.
#' 
#' @param x data.frame returned by sb2hits
#' @param nm integer number of methods in comparison
#'
#' @return
#' tibble with one (binary) column per method, and a `freq` column, with
#' each row corresponding to a single overlap of methods - methods
#' contained in the overlap are set to 1, those not in the overlap
#' are set to 0 - with the `freq` column containing the number of test
#' statistics in the overlap.
#' 
#' @import dplyr magrittr
#' @importFrom UpSetR Counter
#' @author Patrick Kimes
hits2freq <- function(x, nm) {
    UpSetR:::Counter(x, nm, 1, names(x), nintersections = 2^nm,
                     mbar_color = "gray23",
                     order_mat = "degree", aggregate = "degree",
                     cut = NULL, empty_intersects = TRUE,
                     decrease = TRUE) %>%
        as_tibble() %>%
        select(-x, -color)
}


#' Number of Methods w/ Rejections
#'
#' Helper function to return the number of methods with rejections at
#' a particular alpha level (this helps us determine whether or not to plot the
#' aggregated upset plot - if there aren't at least 2 methods it will throw an
#' error, which is a problem for the null simulations).
#' 
#' @param res standardized metric data.table generated using
#'        standardize_results.
#' @param alpha alpha cutoff
#' @param filterSet which methods to exclude from consideration
#'
#' @author Keegan Korthauer
numberMethodsReject <- function(res, alphacutoff, filterSet) {
    res <- res %>% 
        filter(is.na(param.alpha) | (param.alpha == alphacutoff)) %>%
        filter(!(blabel %in% filterSet)) %>%
        filter(alpha == alphacutoff) %>%
        filter(performanceMetric == "rejections") %>%
        select(blabel, performanceMetric, value) %>%
        group_by(blabel) %>%
        summarize(mean_value = round(mean(value))) %>%
        filter(mean_value > 0)
    return(nrow(res))
}
