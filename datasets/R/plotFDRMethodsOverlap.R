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

