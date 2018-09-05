#' Histograms stratified by covariate value quantiles
#' 
#' Takes in a data frame and a covariate of interest, and outputs a cowplot
#' object that displays multiple histograms of p-values: one overall, and 
#' one for each quantile bin (default 3 bins). Also takes as input the bin 
#' width and the max y value for the density.
#'   
#' @param dat data.frame with one row for each hypothesis test, and includes
#'  a column of p-values, and a column of covariate values.
#' @param pvalue character that defines the name of the column in `dat` that 
#'  contains the p-values
#' @param covariate character that defines the name of the column in `dat` that 
#'  contains the covariate
#' @param binwidth a numeric value passed as the parameter for geom_hist    
#' @param maxy a numeric value passed as the upper `limits` parameter to 
#'  `scale_y_continuous`. In other words, this sets the maximum y-value for the 
#'  density histograms (so that all histograms are on the same scale)
#' @param numQ an integer value that specifies the number of quantile bins to 
#'  split on the covariate. Will result in `numQ` + 1 histograms generated. This
#'  must be less than 6.
#' @param main a character string that will be used as the main plot title
#' @param axislabs a logical whether to include x and y axis labels on each panel
#'        
#' @return a cowplot object 
#' 
#' 
#' @author Keegan Korthauer     
strat_hist <- function(dat, pvalue, covariate, binwidth=0.025, maxy=3, numQ=3,
                       main = "", axislabs = TRUE){

  library(cowplot)
  library(ggplot2)
  library(magrittr)
  library(dplyr)
  
  # check numQ input
  if (numQ > 5){
    stop("Please specify a valid value of numQ")
  }
  
  # check for pvalue and covariate cols
  if (!pvalue %in% colnames(dat) | !covariate %in% colnames(dat)){
    stop("pvalue and covariate must be variables that define column names in dat")
  }
  
  plotOne <- function(d, pvalue, covariate, title="", axislabs){
    p<- ggplot(d, aes(x=get(pvalue))) + 
      geom_histogram(binwidth = binwidth, boundary = 0, 
                   colour="grey", fill="lightgrey") +
      aes(y=..density..)+
      theme_classic() +
      theme(plot.subtitle = element_text(size = 11, margin = margin(b = -10), hjust = 0.5),
            plot.title = element_text(size = 11),
            axis.title = element_text(size = 11)) +
      xlim(0,1) +
      scale_y_continuous(expand = c(0.02, 0)) +
      coord_cartesian(ylim=c(0, maxy)) + 
      xlab("") +
      ylab("") +
      labs(title = "      ", subtitle = title) 
    
    if (axislabs)
      p <- p + xlab("p-value") + ylab("Density")
    
    p
  }
  
  gglist <- vector("list", numQ+1)
  
  for (q in 1:(numQ+1)){
    if(q==(numQ+1)){
      gglist[[q]] <- plotOne( dat, pvalue=pvalue, covariate=covariate, title="All",
                              axislabs = axislabs)
    }else{
      dat.strat <- dat %>% 
        dplyr::filter(rank(get(covariate), ties="first") <= quantile(rank(get(covariate), ties="first"), q/numQ) &
               rank(get(covariate), ties="first") >= quantile(rank(get(covariate), ties="first"), (q-1)/numQ) )
      gglist[[q]] <- plotOne(dat.strat, pvalue=pvalue, covariate=covariate, title=paste0("Group ", q),
                             axislabs = axislabs)
    }
  }
   
  if(numQ <= 3){
    numrows <- 1
  }else{
    numrows <- 2
  }
  
  # put last first
  gglist <- gglist[c(numQ+1, 1:numQ)]
  
  gg_stratified <- plot_grid(plotlist=gglist,
                             nrow=numrows)
  
  gg_stratified <- plot_grid(ggplot(), gg_stratified,
                             rel_heights = c(0.075,1),
                             nrow=2) +
                     draw_label(main, 0.5, 0.96)
  return(gg_stratified)
}
