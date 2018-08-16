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
#' @param supplementary Logical. Is the figure a supplementary figure?
#'        
#' @return a ggplot2 object
#' 
#' 
#' @author Alejandro Reyes
rejection_scatter_bins <- function( sb, covariate, threshold=NULL, bins= 4, ncol_facet=2, 
                                    as_fraction=FALSE, supplementary=TRUE,
                                    palette = candycols ){
  stopifnot(is(sb, "SummarizedExperiment"))
  alphas <- as.numeric( as.character( unique( colData(sb)$param.alpha ) ) )
  alphas <- alphas[!is.na(alphas)]
  if( !is.null(threshold) ){
    sb <- sb[rowData(sb)[[covariate]] > threshold,]
  }
  if(sum(colSums(is.na(assays(sb)$qvalue)) == nrow(sb)) > 0){
    miss <- which(colSums(is.na(assays(sb)$qvalue)) == nrow(sb)) 
    sb <- sb[,-miss]
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
                                    function(x){paste( signif( range(x), 2 ), collapse=" - ")} )
  plotDF <- dataPerBin %>% dplyr:::filter( !(grepl("ihw", blabel) & param.alpha != alpha ) ) %>%
    dplyr:::mutate( blabel=gsub("(ihw)-a\\d+", "\\1", blabel ) ) %>%
#    dplyr:::select( blabel, key, value, assay, performanceMetric, alpha, bin ) %>%
    dplyr:::filter( blabel!="unadjusted", performanceMetric == "rejections" ) 
  if( !supplementary ){
    plotDF <- plotDF %>%
      dplyr:::mutate( param.smooth.df=gsub("L", "", param.smooth.df ) ) %>%
      dplyr:::filter( !( grepl("bl-df", blabel) & 
                           as.numeric(as.character(param.smooth.df) != 3 ) ) ) %>%
      dplyr:::mutate( Method=gsub("-df03", "", blabel))
  }
  
  # add color palette
  plotDF <- dplyr::left_join(plotDF, palette, by="Method") 
  
  col <- as.character(plotDF$col)
  names(col) <- as.character(plotDF$Method)
  
  lty <- as.character(plotDF$lty)
  names(lty) <- as.character(plotDF$Method)

  plotDF %>%
    ggplot( aes(alpha, value, col=Method) ) +
    geom_line(alpha = 3/4, aes(linetype=Method)) +  
    geom_point(alpha = 3/4, show.legend = FALSE) + 
    facet_wrap(~bin, ncol=ncol_facet ) +
    xlab(expression(paste("Nominal"~alpha))) +
    scale_color_manual(values=col) +
    scale_linetype_manual(values = lty) +
    ylab(yl) 
}
