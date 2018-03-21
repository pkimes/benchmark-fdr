#' Boxplots of covariates for rejections unique to methods subsets
#'   
#' @param object A SummarizedBenchmark object
#' @param alpha An alpha value.
#' @param nsets Number of method sets to plot. Method sets are ordered decreasingly according 
#' to the number of rejections specific to the sets of methods.
#' @param methods A character vector of methods to consider (must be a level of colData(object)$blabel).
#' @param trans A function to transform the y-axis of the boxplot
#' @param maxNum An integer that represents the maximum number of methods per set
#'        
#' @return an ggplot2 boxplot
#' 
#' @author Alejandro Reyes
plotCovariateBoxplots <- function( object, alpha, nsets=8, methods=NULL, 
                                   maxNum=length(methods), 
                                   trans=function(x){x} ){
  stopifnot( is(object, "SummarizedBenchmark") )
  if( is.null( methods ) ){
    methods <- colData(object)$blabel
  }
  mat <- ( assays( object )[["qvalue"]][,methods] < alpha )
  covariate <- mcols( object )$ind_covariate
  allLengths <- lapply( seq_len( maxNum ), function(n){
    combs <- combn( methods, n )
    indCov <- lapply( seq_len( ncol( combs ) ), function(x){
      detectedInComb <- rowSums( mat[,combs[,x], drop=FALSE], na.rm=TRUE ) == n & 
        !(rowSums( mat[,methods,drop=FALSE], na.rm=TRUE ) > n)
      covariate[detectedInComb]
    } )
    combs <- apply( combs, 2, function(x){ paste( gsub("-", "", x), collapse="-" ) } )
    names(indCov) <- combs
    indCov[sapply( indCov, length ) > 0]
  } )
  allLengths <- unlist( allLengths, recursive=FALSE )
  allLenghts <- allLengths[order( sapply( allLengths, length ), decreasing = TRUE)]
  allLenghts <- head( allLenghts, nsets )
  sizesAll <- unlist( allLenghts )
  dfCov <- data.frame(
    covariate = sizesAll,
    methods=factor( rep( names( allLenghts ), sapply( allLenghts, length ) ), 
                    levels=names(allLenghts)  ) )
  pl <- ggplot( dfCov, aes( methods, trans(covariate), col=methods) ) +
    geom_boxplot( ) + theme(legend.position="top", 
                            axis.text.x=element_blank(), 
                            legend.direction="vertical") +
    labs(col="", y="Covariate", x="") + guides(color=guide_legend(ncol=2))
  pl
}
