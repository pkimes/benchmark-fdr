#' Differential Abundance Testing in Microbiome Data
#'
#' This wrapper function performs differential testing of
#' OTU (relative) abundances between case and control samples
#' using the Wilcoxon rank sum test (equiv. Mann-Whitney test)
#' for each OTU. Results are reported as a data.frame with
#' one row per OTU.
#' 
#' @param abundance matrix of relative abundances with rows and
#'        columns corresponding to OTUs and samples.
#' @param shift psuedocount that was added to all relative abundance
#'        values in \code{abundance} to handle zero counts.
#' @param is_case logical vector of length equal to the number of
#'        columns in \code{abundance} specifying whether each sample
#'        is a case (TRUE) or control (FALSE).
#'
#' @return
#' data.frame of differential testing results.
#'
#' @author Claire Duvallet
test_microbiome <- function(abundance, shift, is_case) {

    casedf <- abundance[, is_case]
    ctrldf <- abundance[, !is_case]

    ## Calculate pvalues, effects, and stderr
    res <- data.frame(otu = rownames(abundance), 
                      pval = NA, wilcox_teststat = NA,
                      SE = NA, effect_size = NA, 
                      mean_abun = NA, mean_abun_present = NA,
                      ubiquity = NA)
    
    for (i in seq_len(nrow(abundance))) {
        ## wilcoxon p value
        w <- wilcox.test(casedf[i, ], ctrldf[i, ])
        res$pval[i] <- w$p.value
        res$wilcox_teststat[i] <- w$statistic
                
        ## standard error of the OTU abundance, across all samples
        res$SE[i] <- sd(abundance[i, ]) / sqrt(length(abundance[i, ]))
        
        ## mean OTU abundance across all samples (after removing pseudo-count)
        res$mean_abun[i] <- mean(abundance[i, ]) - shift
        
        ## mean OTU abundance across only samples with the OTU present
        res$mean_abun_present[i] <-
            sum(abundance[i, ] - shift) / sum(abundance[i, ] > shift)
        
        ## ubiquity of OTU across all samples
        res$ubiquity[i] <- sum(abundance[i, ] > shift) / length(abundance[i, ])
        
        ## effect (logfold difference)
        res$effect_size[i] <- log( mean(ctrldf[i, ]) / mean(casedf[i, ]) )
    }

    ## compute test-statistic from p-value for ASH
    res <- dplyr::mutate(res,
                         test_statistic = qnorm(exp(log(pval) - log(2)),
                                                lower.tail=FALSE),
                         test_statistic = test_statistic * sign(effect_size))
}

test_microbiome_corr <- function(abundance, meta_col, shift) {
  
  ## Calculate correlation and pvalues
  res <- data.frame(otu = rownames(abundance), 
                    pval = NA, spearman_teststat = NA,
                    SE = NA, effect_size = NA, 
                    mean_abun = NA, mean_abun_present = NA,
                    ubiquity = NA)
  
  for (i in seq_len(nrow(abundance))) {
    ## spearman p value
    corr <- cor.test(abundance[i, ], meta_col, 
                     method='spearman', exact = FALSE, alternative = "two.sided")
    res$pval[i] <- corr$p.value
    res$spearman_teststat[i] <- corr$statistic
    
    ## standard error of the OTU abundance, across all samples
    res$SE[i] <- sd(abundance[i, ]) / sqrt(length(abundance[i, ]))
    
    ## mean OTU abundance across all samples (after removing pseudo-count)
    res$mean_abun[i] <- mean(abundance[i, ]) - shift
    
    ## mean OTU abundance across only samples with the OTU present
    res$mean_abun_present[i] <-
      sum(abundance[i, ] - shift) / sum(abundance[i, ] > shift)
    
    ## ubiquity of OTU across all samples
    res$ubiquity[i] <- sum(abundance[i, ] > shift) / length(abundance[i, ])
    
  }
  
  ## compute test-statistic from p-value for ASH
  res <- dplyr::mutate(res,
                       test_statistic = qnorm(exp(log(pval) - log(2)),
                                              lower.tail=FALSE),
                       test_statistic = test_statistic * sign(effect_size))
}
  

remove_shallow_smpls <- function(df, n_reads) {
    ## df has OTUs in rows and samples in columns
    return(df[, colSums(df) > n_reads])
}


remove_shallow_otus <- function(df, n_reads) {
    return(df[rowSums(df) > n_reads, ])
}

remove_rare_otus <- function(df, perc_samples){
  return(df[rowSums(df > 0) / dim(df)[2] > perc_samples, ])
}
