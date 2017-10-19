#' Rank scatter plots
#' 
#' Function to plot a scatter plot of the covariate rank by the -log10 p-value
#'
#' The rank scatter function takes as input a data frame and the 
#' covariate of interest
#' and outputs a ggplot object that displays a scatter plot of the covariate
#' rank versus the -log10 p-value.
#' It also takes as input the number of bins for geom_hex
#'
#' @param dat data.frame with one row for each hypothesis test, and includes
#'  a column of p-values, and a column of covariate values.
#' @param pval character that defines the name of the column in `dat` that 
#'  contains the p-values
#' @param covariate character that defines the name of the column in `dat` that 
#'  contains the covariate
#' @param bins a numeric value passed as the bins parameter for geom_hex
#' 
#' @return a ggplot object 
#' 
#' @author Keegan Korthauer      
rank_scatter <- function(dat, pval, covariate, bins=100){
  gg_scat <- ggplot(dat, aes(y=-log10(get("pval")), 
                             x=rank(get(covariate), ties="first")/nrow(dat))) +
    geom_hex(bins = bins) +
    ylab(expression(-log[10]~p)) +
    xlab("Covariate Rank")
  theme_classic()
  return(gg_scat)
}
