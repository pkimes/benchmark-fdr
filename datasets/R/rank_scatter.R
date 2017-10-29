#' Scatter plots (p-value vs covariate)
#' 
#' Function to plot a scatter plot of the covariate by the -log10 p-value
#'
#' This function takes as input a data frame and the covariate of interest
#' and outputs a ggplot object that displays a scatter plot of the covariate
#' rank versus the -log10 p-value.
#' It also takes as input the number of bins for geom_hex
#'
#' @param dat data.frame with one row for each hypothesis test, and includes
#'  a column of p-values, and a column of covariate values.
#' @param pvalue character that defines the name of the column in `dat` that 
#'  contains the p-values
#' @param covariate character that defines the name of the column in `dat` that 
#'  contains the covariate
#' @param bins a numeric value passed as the bins parameter for geom_hex
#' @param funx function to transform the x-axis
#' @param funfill function from the scales package to transform the color scale 
#' 
#' @return a ggplot object 
#' 
#' @author Keegan Korthauer      
rank_scatter <- function( dat, pvalue, covariate,bins=100,
                         funx=function(x){rank(x, ties="first")/nrow(dat)},
                         funfill=NULL )
{
  gg_scat <- ggplot(dat, aes(y=-log10(get(pvalue)), x=funx(get(covariate)))) +
    geom_hex(bins = bins)+
    ylab(expression(-log[10]~p)) +
    xlab("Covariate") +
    theme_classic()
  if( !is.null( funfill ) ){
    gg_scat <- gg_scat + scale_fill_continuous(trans=funfill)
  }
  return(gg_scat)
}
