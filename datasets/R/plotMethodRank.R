#' Function to create a heatmap of method rankings for different case studies,
#' datasets, or covariates. 
#' 
#' Takes as input a vector of filepaths, each one pointing to a SummarizedBenchmark
#' results object. Output is a ggplot object with the heatmap. Each file 
#' specified will represent one column of the heatmap (rows are methods).
#' 
#' @param objects a character vector containing the full filepaths to 
#'  SummarizedBenchmark results objects.
#' @param colLabels a character vector containing the labels to use for the 
#'  heatmap columns (same length and order as `objects`). 
#' @param alpha numeric value indicating which alpha value to use for IHW 
#'  (needs to be one of the values specified in the bench object). Default 0.10.
#' @param colorLow character representing the color to be used for the lowest 
#'  value
#' @param colorHigh character representing the color to be used for the highest
#'  value
#' @param colorNA character representing the color to be used for NA (i.e. 
#'  for methods that weren't used in the particular benchmark)
#' @param linePlot logical indicating whether to plot a line plot in lieu of 
#'  a heatmap  
#' @param excludeMethods character vector containing the names of methods to
#'  exclude from the heatmap. Default is to exclude Scott methods which are 
#'  only present in GWAS.
#' @param propMaxRejections logical indicating whether the proportion of max
#'  number of rejections is plotted (instead of the rank). Default is TRUE. 
#' @param xlab character with the x label (default is "Case Study")   

plotMethodRanks <- function(objects, colLabels, alpha = 0.10, 
                            colorLow = "navy", colorHigh = "yellow",
                            colorNA = "white", xlab = "Case Study",
                            propMaxRejections = TRUE,
                            linePlot = FALSE,
                            excludeMethods = c("scott-theoretical",
                                               "scott-empirical")){

  # create tidy data frame where each row is a method / dataset observation
  # of a rank 
  ranks <- data.frame()
  
  for (i in seq_along(objects)){
    x <- readRDS(objects[i])
    assayNames(x) <- "qvalue"
    x <- addDefaultMetrics( x )
    
    hasResults <- apply(!is.na( assays(x)[["qvalue"]] ), 2, sum)
    
    if (sum(hasResults) > 0){
      tmp <- estimatePerformanceMetrics(x, alpha, tidy=TRUE) %>%
        filter( performanceMetric == "rejections") %>%
        dplyr::rename( method = blabel) %>%
        filter( is.na(param.alpha) | (param.alpha == alpha)) %>%
        filter( is.na(param.smooth.df) | (param.smooth.df == "3L")) %>%
        filter( !method == "unadjusted") %>%
        select( method, value ) %>%
        dplyr::rename( nrejects = value) %>%
        mutate( method = gsub("-df03", "", method)) %>%
        mutate( method = gsub("(-a)(.*)", "", method)) %>%
        na.omit() %>%
        mutate( rank = rank(nrejects) / n(),
                propMaxRej = nrejects / max(nrejects, na.rm = TRUE)) %>%
        mutate(casestudy = colLabels[i])
      
      ranks <- rbind(ranks, tmp)
    }
  }
  
  # exclude methods 
  if(!is.null(excludeMethods)){
    ranks <- ranks %>% filter( !(method %in% excludeMethods))
  }
  
  # add NAs for missing combinations of methods and casestudies
  ranks <- ranks %>% complete(method, casestudy)
  
  # average over datasets within case study
  ranks_avg <- ranks %>% 
    group_by( casestudy, method) %>%
    summarize( mean_rank = mean(rank, na.rm = TRUE),
               mean_prop = mean(propMaxRej, na.rm = TRUE)) 
  
  case_dat <- ranks_avg %>%
    group_by( casestudy) %>%
    summarize( nmethods = sum(!is.na(mean_rank))) 
  
  method_dat <- ranks_avg %>%
    group_by( method) %>%
    summarize( meta_rank = mean(mean_rank, na.rm = TRUE),
               meta_prop = mean(mean_prop, na.rm = TRUE))
  
  ranks_avg <- left_join(ranks_avg, method_dat, by ="method")
  ranks_avg <- left_join(ranks_avg, case_dat, by ="casestudy")
  
  # reorder columns
  ranks_avg$casestudy <- factor(ranks_avg$casestudy, 
                                levels=unique(ranks_avg$casestudy[order(ranks_avg$nmethods,
                                                                        decreasing= TRUE)]))
  
  if (!linePlot){
    if (propMaxRejections){
      # reorder rows
      ranks_avg$method <- factor(ranks_avg$method,
                                 levels=unique(ranks_avg$method[order(ranks_avg$meta_prop)]))
      
      # heatmap : rows method, cols casestudy
      Fig <- ggplot(ranks_avg, aes(x = casestudy, y = method, fill = mean_prop)) + 
        geom_raster() + 
        scale_fill_gradient(low = colorLow, high = colorHigh, na.value = colorNA) +
        theme_bw() +
        xlab(xlab) +
        ylab("Method") +
        labs(fill = "Mean proportion\nMax rejections")
    }else{
      # reorder rows
      ranks_avg$method <- factor(ranks_avg$method,
                                 levels=unique(ranks_avg$method[order(ranks_avg$meta_rank)]))
      
      Fig <- ggplot(ranks_avg, aes(x = casestudy, y = method, fill = mean_rank)) + 
        geom_raster() + 
        scale_fill_gradient(low = colorLow, high = colorHigh, na.value = colorNA) +
        theme_bw() +
        xlab(xlab) +
        ylab("Method") +
        labs(fill = "Mean rank")
    }
  }else{
    # reorder rows
    ranks_avg$method <- factor(ranks_avg$method,
                               levels=rev(unique(ranks_avg$method[order(ranks_avg$meta_rank)])))
    Fig <- ggplot( ranks_avg, aes(x=method, y=mean_rank, 
                                    group=casestudy, col=casestudy) ) +
      geom_line() + geom_point() +
      labs(y="Mean rank", x="", col=xlab) +
      theme(axis.text.x=element_text(angle=25, vjust=1, hjust=1))
  }  

  return(Fig)
}
