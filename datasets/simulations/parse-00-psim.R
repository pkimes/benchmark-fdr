#' Compute All Performance Metrics for Simulations
#'
#' @author Patrick Kimes

## Workspace Setup
library("SummarizedBenchmark")
library("tidyverse")
library("magrittr")

## helper function to loop over 
source("../R/plotsim_standardize.R")

## sequence of alpha values to consider
alpha_set <- c(seq(0.01, 0.10, 0.01), 0.15, 0.20, 0.25)

## look at each individual simulation separately
df <- list.files("data-psim/M100", "*-M100\\.rds", recursive = TRUE, full.names = TRUE)

## parse each setting-block separately and save metrics in single file
for (f in df) {
    outf <- gsub("data-psim", "data-psim-parsed", f)
    if (!file.exists(outf)) {
        print(paste0("running: ", basename(f)))
        dir.create(dirname(outf))
        ## read process data
        res <- readRDS(f)
        tsb <- plotsim_standardize(res, alpha = alpha_set) %>%
            select(blabel, performanceMetric, value, rep, alpha)
        saveRDS(tsb, outf)
        print(paste0("done with: ", basename(f))) 
    } else {
        print(paste0("skipping: ", basename(f)))
    }
}


