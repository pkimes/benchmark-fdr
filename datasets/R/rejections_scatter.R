#' Scatterplot of nominal alpha vs rejections
#' 
#' This function inputs a SummarizedBenchmark object build by the
#' common benchDesign in this repository. It calculates the number 
#' of rejections for several nominal alpha values and plots them. 
#'   
#' @param sb A SummarizedBenchmark object
#' @param asFraction Logical. Whether to plot the fraction of hypotheses
#' that were rejected (TRUE) or the absolute number of rejections (FALSE).
#'        
#' @return a ggplot2 object
#' 
#' 
#' @author Alejandro Reyes
rejections_scatter <- function( sb, asFraction=FALSE ){
  stopifnot( is(sb, "SummarizedBenchmark" ) )
  alphas <- unique( as.numeric( as.character( colData( sb )$param.alpha ) ) )
  alphas <- alphas[!is.na(alphas)]
  if( asFraction ){
    deno <- nrow(sb)
    yl <- "Fraction of hypotheses rejected"
  }else{
    deno <- 1
    yl <- "Number of rejections"
  }
  estimatePerformanceMetrics( sb, alpha=alphas, tidy=TRUE ) %>%
    dplyr:::filter( !(pkg_name == "IHW" & param.alpha != alpha ) ) %>%
    dplyr:::mutate( blabel=gsub("(ihw)-a\\d+", "\\1", blabel ) ) %>%
    dplyr:::select( blabel, key, value, assay, performanceMetric, alpha ) %>%
    dplyr:::filter( blabel!="unadjusted") %>%
    ggplot( aes(alpha, value/deno, col=blabel) ) +
    geom_line() + geom_point() +
    xlab(expression(paste("Nominal"~alpha))) +
    ylab(yl)
}