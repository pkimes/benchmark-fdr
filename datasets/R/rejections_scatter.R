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
rejections_scatter <- function( sb, as_fraction=FALSE ){
  stopifnot( is(sb, "SummarizedBenchmark" ) )
  alphas <- unique( as.numeric( as.character( colData( sb )$param.alpha ) ) )
  alphas <- alphas[!is.na(alphas)]
  if( as_fraction ){
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
    dplyr:::filter( blabel!="unadjusted", performanceMetric == "rejections" ) %>%
    ggplot( aes(alpha, value/deno, col=blabel) ) +
    geom_line() + geom_point() +
    xlab(expression(paste("Nominal"~alpha))) +
    ylab(yl)
}

#' Scatterplot of nominal alpha vs rejections stratified by covariate bins
#' 
#' Similar to the function rejections_scatter, but it stratifies the scatter 
#' plot by covariate bins.  
#'   
#' @param sb A SummarizedBenchmark object
#' @param as_fraction Logical. Whether to plot the fraction of hypotheses
#' that were rejected (TRUE) or the absolute number of rejections (FALSE).
#' @param threshold A numeric threshold. If specified, only test where the covariate
#' is larger than this number are considered. 
#' @param bins Number of bins.
#' @param ncol_facet Number of columns for facet_grid
#' @param covariate Character object that represents the name of the rowData
#' column containing the independent covariate to bin on
#'        
#' @return a ggplot2 object
#' 
#' 
#' @author Alejandro Reyes
rejection_scatter_bins <- function( sb, covariate, threshold=NULL, bins= 4, ncol_facet=2, as_fraction=FALSE ){
  stopifnot(is(sb, "SummarizedExperiment"))
  alphas <- as.numeric( as.character( unique( colData(sb)$param.alpha ) ) )
  alphas <- alphas[!is.na(alphas)]
  if( !is.null(threshold) ){
    sb <- sb[rowData(sb)[[covariate]] > threshold,]
  }
  rowData(sb)$bin <- cut(rank( rowData(sb)[[covariate]], ties="first"),
                         quantile( rank(rowData(sb)[[covariate]], ties="first"), seq(0, 1, length.out=bins+1) ),
                         include.lowest=TRUE)
  dataPerBin <- lapply( levels( rowData(sb)$bin), function(x){
    inBin <- rowData(sb)$bin == x
    res <- estimatePerformanceMetrics( sb[inBin,], tidy=TRUE, alpha=alphas )
    if( as_fraction ){
      res <- res %>% mutate( value=value/sum(inBin) )
    }
    cbind( res, bin=x )
  } )
  yl <- ifelse(as_fraction, "Fraction of hypotheses rejected", "Number of rejections" )
  dataPerBin <- do.call(rbind, dataPerBin)
  levels(dataPerBin$bin) <- tapply( rowData(sb)[[covariate]], 
                                    rowData(sb)$bin, 
                                    function(x){paste(range(x), collapse=" - ")} )
  dataPerBin %>% dplyr:::filter( !(pkg_name == "IHW" & param.alpha != alpha ) ) %>%
    dplyr:::mutate( blabel=gsub("(ihw)-a\\d+", "\\1", blabel ) ) %>%
    dplyr:::select( blabel, key, value, assay, performanceMetric, alpha, bin ) %>%
    dplyr:::filter( blabel!="unadjusted", performanceMetric == "rejections" ) %>%
    ggplot( aes(alpha, value, col=blabel) ) +
    geom_line() + geom_point() + facet_wrap(~bin, ncol=ncol_facet ) +
    xlab(expression(paste("Nominal"~alpha))) +
    ylab(yl)
}
