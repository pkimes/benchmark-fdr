#' Scatterplot of nominal alpha vs rejections
#' 
#' This function inputs a SummarizedBenchmark object build by the
#' common benchDesign in this repository. It calculates the number 
#' of rejections for several nominal alpha values and plots them. 
#'   
#' @param sb A SummarizedBenchmark object
#' @param asFraction Logical. Whether to plot the fraction of hypotheses
#' that were rejected (TRUE) or the absolute number of rejections (FALSE).
#' @param supplementary Logical. Is the figure a supplementary figure?
#'        
#' @return a ggplot2 object
#' 
#' 
#' @author Alejandro Reyes
rejections_scatter <- function( sb, as_fraction=FALSE, supplementary=TRUE,
                                palette = candycols){
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
  if(sum(colSums(is.na(assays(sb)$qvalue)) == nrow(sb)) > 0){
    miss <- which(colSums(is.na(assays(sb)$qvalue)) == nrow(sb)) 
    sb <- sb[,-miss]
  }
  
  plotDF <- estimatePerformanceMetrics( sb, alpha=alphas, tidy=TRUE ) %>%
    dplyr:::filter( !(grepl("ihw", blabel) & param.alpha != alpha ) ) %>%
    dplyr:::mutate( blabel=gsub("(ihw)-a\\d+", "\\1", blabel ) ) %>%
    #dplyr:::select( blabel, key, value, assay, performanceMetric, alpha ) %>%
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
    ggplot( aes(alpha, value/deno, col=Method) ) +
    geom_line(alpha = 3/4, aes(linetype=Method)) + 
    geom_point(alpha = 3/4, show.legend = FALSE) +
    xlab(expression(paste("Nominal"~alpha))) +
    scale_color_manual(values=col) +
    scale_linetype_manual(values = lty) +
    ylab(yl)
}
