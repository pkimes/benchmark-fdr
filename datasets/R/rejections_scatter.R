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
rejections_scatter <- function( sb, as_fraction=FALSE, supplementary=TRUE ){
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
      dplyr:::mutate( blabel=gsub("-df03", "", blabel))
  }
  plotDF %>%
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
#' @param supplementary Logical. Is the figure a supplementary figure?
#'        
#' @return a ggplot2 object
#' 
#' 
#' @author Alejandro Reyes
rejection_scatter_bins <- function( sb, covariate, threshold=NULL, bins= 4, ncol_facet=2, 
                                    as_fraction=FALSE, supplementary=TRUE ){
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
                                    function(x){paste( round( range(x), 2 ), collapse=" - ")} )
  plotDF <- dataPerBin %>% dplyr:::filter( !(grepl("ihw", blabel) & param.alpha != alpha ) ) %>%
    dplyr:::mutate( blabel=gsub("(ihw)-a\\d+", "\\1", blabel ) ) %>%
#    dplyr:::select( blabel, key, value, assay, performanceMetric, alpha, bin ) %>%
    dplyr:::filter( blabel!="unadjusted", performanceMetric == "rejections" ) 
  if( !supplementary ){
    plotDF <- plotDF %>%
      dplyr:::mutate( param.smooth.df=gsub("L", "", param.smooth.df ) ) %>%
      dplyr:::filter( !( grepl("bl-df", blabel) & 
                           as.numeric(as.character(param.smooth.df) != 3 ) ) ) %>%
      dplyr:::mutate( blabel=gsub("-df03", "", blabel))
  }
  plotDF %>%
    ggplot( aes(alpha, value, col=blabel) ) +
    geom_line() + geom_point() + facet_wrap(~bin, ncol=ncol_facet ) +
    xlab(expression(paste("Nominal"~alpha))) +
    ylab(yl)
}

#' Plot overlaps between FDR methods
#'   
#' @param object A SummarizedBenchmark object
#' @param alpha An alpha value.
#' @param supplementary Logical. Is the figure a supplementary figure?
#'        
#' @return an upsetr plot
#' 
#' @author Alejandro Reyes
plotFDRMethodsOverlap <- function( object, supplementary=TRUE, alpha=0.1, ... ){
  stopifnot( is( object, "SummarizedBenchmark" ) )
  stopifnot( any( colData(object)$param.alpha == alpha, na.rm=TRUE) )
  object <- object[,!( grepl("^ihw", as.character( colData( object )$blabel ) ) & colData( object )$param.alpha != alpha )]
  colData(object)$blabel <- gsub("(ihw)-.*", "\\1", colData( object )$blabel)
  qvals <- assays( object )[["qvalue"]]
  object <- object[,!apply( is.na( assays( object )[["qvalue"]] ), 2, all )]
  object <- object[,colData(object)$blabel != "unadjusted"]
  if( !supplementary ){
    object <- object[,!( grepl("bl", colData( object )$blabel) & as.numeric( gsub( "L", "", colData( object )$param.smooth.df ) ) != 3 )]
    colData(object)$blabel <- gsub( "(bl)-.*", "\\1", colData(object)$blabel )
  }
  colnames(object) <- colData( object )$blabel
  plotMethodsOverlap( object, alpha=alpha, ... )
}


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
