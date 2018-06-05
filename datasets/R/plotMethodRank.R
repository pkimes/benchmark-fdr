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
#' @param fill character indicating the outcome variable of interest. Default
#'  is 'propMaxRejections'. For 'FDR' and 'TPR' sb objects need to contain 
#'  these performance metrics.
#' @param xlab character with the x label (default is "Case Study")   
#' @param tableOnly logical whether or not to return a tibble of the summarized
#'  results instead of a plot.
#' @param rowOrder numerical vector with order of rows/methods (default is 
#'  NULL, which means rows will be ordered by mean fill value).
#' @param Nlabel logical whether or not to include the number of datasets 
#'  averaged over in the x-axis casestudy label. Default is TRUE. If FALSE
#'  assumes each column represents one benchmark object (legend labels will not 
#'  include the word "Mean").
#' @param annotate character indicating what should be plotted as text labels
#'  (heatmap plot only). Default is proportion of total possible rejections 
#'  (propPossible). Could also be another variable (from \code{fill} choices) 
#'  or NULL (for no text label annotations).
#' @author Keegan Korthauer  

plotMethodRanks <- function(objects, colLabels, alpha = 0.10, 
                            colorLow = "navy", colorHigh = "yellow",
                            colorNA = "white", xlab = "Case Study",
                            fill = c("propMaxRejections", "meanRank", 
                                     "FDR", "TPR"),
                            linePlot = FALSE,
                            excludeMethods = c("scott-theoretical",
                                               "scott-empirical"),
                            rowOrder = NULL,
                            tableOnly = FALSE,
                            Nlabel = TRUE,
                            annotate = "propPossible"){
  fill <- match.arg(fill)
  if (is.list(readRDS(objects[1]))){
    ranks <- data.frame()
    
    for (l in seq_along(objects)){
      x <- readRDS(objects[l])
      tmp <- tidy_df(x, colLabels = rep(colLabels[l], length(x)), 
                     fill, annotate) %>%
        group_by(method, casestudy) %>%
        summarize(nrejects = mean(nrejects),
                  rank = mean(rank),
                  propMaxRej = mean(propMaxRej),
                  propPossible= mean(propPossible))
      ranks <- rbind(ranks, as.data.frame(tmp))
    }
    
  }else{
    ranks <- tidy_df(objects, colLabels, fill, annotate)
  }
  
  # exclude methods 
  if(!is.null(excludeMethods)){
    ranks <- ranks %>% dplyr::filter( !(method %in% excludeMethods))
  }
  
  # add NAs for missing combinations of methods and casestudies
  ranks <- ranks %>% complete(method, casestudy)
  
  # average over datasets within case study
  if (fill %in% c("FDR", "TPR")){
    annot <- "nrejects"
  }else{
    annot <- annotate
  }
  
  if(is.null(annot))
    annot <- "propPossible"
  ranks_avg <- ranks %>% 
    group_by( casestudy, method) %>%
    summarize( mean_rank = mean(rank, na.rm = TRUE),
               mean_prop = mean(propMaxRej, na.rm = TRUE),
               mean_nrej = mean(nrejects, na.rm = TRUE),
               min_prop = min(propMaxRej, na.rm = TRUE),
               max_prop = max(propMaxRej, na.rm = TRUE),
               topLayer = mean(get(annot), na.rm = TRUE),
               nsets = sum(!is.na(nrejects)))
  
  if(Nlabel){
    if(!is.list(readRDS(objects[1]))){
      ranks_avg <- ranks_avg %>% 
        group_by( casestudy ) %>%
        mutate( studyname = paste0(casestudy, "(", max(nsets), ")"),
                nsets = ifelse(nsets==0, NA, nsets))
    }else{
      ranks_avg <- ranks_avg %>% 
        group_by( casestudy ) %>%
        mutate( studyname = paste0(casestudy, "(", 
                                 length(readRDS(objects[1])), ")"),
              nsets = ifelse(nsets==0, NA, nsets))
    } 
  }else{
    ranks_avg <- ranks_avg %>% mutate(studyname = casestudy)
  }
  
  case_dat <- ranks_avg %>%
    group_by( casestudy) %>%
    summarize( nmethods = sum(!is.na(mean_rank))) 
  
  method_dat <- ranks_avg %>%
    group_by( method) %>%
    summarize( meta_rank = mean(mean_rank, na.rm = TRUE),
               meta_prop = mean(mean_prop, na.rm = TRUE),
               meta_nrej = mean(mean_nrej, na.rm = TRUE))
  
  ranks_avg <- left_join(ranks_avg, method_dat, by ="method")
  ranks_avg <- left_join(ranks_avg, case_dat, by ="casestudy")
  
  # reorder columns
  ranks_avg$casestudy <- factor(ranks_avg$casestudy, 
                                levels=unique(ranks_avg$casestudy[order(ranks_avg$nmethods,
                                                                        decreasing= TRUE)]))
  ranks_avg$studyname <- factor(ranks_avg$studyname, 
                                levels=unique(ranks_avg$studyname[order(ranks_avg$nmethods,
                                                                        decreasing= TRUE)]))
  
  if(!is.null(rowOrder))
    ranks_avg$method <- factor(ranks_avg$method,
                               levels=rowOrder)
  if (tableOnly){
    ranks_avg$method <- factor(ranks_avg$method,
                               levels=unique(ranks_avg$method[order(ranks_avg$meta_prop)]))
    return(ranks_avg)
  }else if (!linePlot){
    if (fill == "propMaxRejections"){
      # reorder rows
      if (is.null(rowOrder))
        ranks_avg$method <- factor(ranks_avg$method,
                                   levels=unique(ranks_avg$method[order(ranks_avg$meta_prop)]))
      
      # heatmap : rows method, cols casestudy
      Fig <- ggplot(ranks_avg, aes(x = studyname, y = method, fill = mean_prop)) + 
        geom_raster() + 
        scale_fill_gradient(low = colorLow, high = colorHigh, na.value = colorNA,
                            limits=c(0,1)) +
        theme_bw() +
        xlab(xlab) +
        ylab("Method") +
        labs(fill = "Mean proportion\nmax rejections")
      
      if(!Nlabel)
        Fig <- Fig +
          labs(fill = "Proportion max\nrejections") 
    }else if (fill == "meanRank"){
      # reorder rows
      if (is.null(rowOrder))
        ranks_avg$method <- factor(ranks_avg$method,
                                 levels=unique(ranks_avg$method[order(ranks_avg$meta_rank)]))
      
      Fig <- ggplot(ranks_avg, aes(x = studyname, y = method, fill = mean_rank)) + 
        geom_raster() + 
        scale_fill_gradient(low = colorLow, high = colorHigh, na.value = colorNA) +
        theme_bw() +
        xlab(xlab) +
        ylab("Method") +
        labs(fill = "Mean rank")
      
      if(!Nlabel)
        Fig <- Fig +
        labs(fill = "Rank") 
    }else{
      # reorder rows
      if (is.null(rowOrder))
        ranks_avg$method <- factor(ranks_avg$method,
                                 levels=unique(ranks_avg$method[order(ranks_avg$meta_nrej)]))
      
      Fig <- ggplot(ranks_avg, aes(x = studyname, y = method, fill = mean_nrej)) + 
        geom_raster() + 
        scale_fill_gradient(low = colorLow, high = colorHigh, na.value = colorNA) +
        theme_bw() +
        xlab(xlab) +
        ylab("Method") +
        labs(fill = paste0("Mean ", fill))
    }
    
    if (!is.null(annotate)){
      Fig <- Fig + 
        geom_text(data = ranks_avg, 
                  aes(label = ifelse(is.na(topLayer), "", 
                                     sprintf("%.3f", round(topLayer,3)))))
    }
    
    return(Fig)
  }else{
    message("Line plot only supported for mean rank")
    # reorder rows
    if (is.null(orderRows))
      ranks_avg$method <- factor(ranks_avg$method,
                               levels=rev(unique(ranks_avg$method[order(ranks_avg$meta_rank)])))
    Fig <- ggplot( ranks_avg, aes(x=method, y=mean_rank, 
                                    group=casestudy, col=casestudy) ) +
      geom_line() + geom_point() +
      labs(y="Mean rank", x="", col=xlab) +
      theme(axis.text.x=element_text(angle=25, vjust=1, hjust=1))
    return(Fig)
  }  
}



tidy_df <- function(objects, colLabels, fill, annotate){
    ## create tidy data frame where each row is a method / dataset observation
    ## of a rank 
    ranks <- data.frame()
    
    for (i in seq_along(objects)){
        if ( class(objects[[i]]) == "character") {
            x <- readRDS(objects[i])
        }else if ( class(objects[[i]]) == "SummarizedBenchmark") {
            x <- objects[[i]]
        }
        
        assayNames(x) <- "qvalue"
        x <- addDefaultMetrics( x )
        
        hasResults <- apply(!is.na( assays(x)[["qvalue"]] ), 2, sum)
        NAmethods <- names(hasResults)[hasResults == 0]
        
        if (fill %in% c("propMaxRejections", "meanRank")){
            pmcol <- "rejections"
        }else{
            pmcol <- fill
        }
        
        if (sum(hasResults) > 0){
            tmp <- estimatePerformanceMetrics(x, alpha, tidy=TRUE)
            if (fill %in% c("TPR", "FDR") && !(fill %in% tmp$performanceMetric))
                stop(fill, " is not found in performanceMetrics")
            
            if (!is.null(annotate)){
                if (annotate %in% c("TPR", "FDR") && !(annotate %in% tmp$performanceMetric))
                    stop(annotate, " is not found in performanceMetrics")
            }
            
            tmp <- tmp %>%
                dplyr::filter( performanceMetric == pmcol) %>%
                dplyr::rename( method = blabel) %>%
                dplyr::filter( is.na(param.alpha) | (param.alpha == alpha)) %>%
                dplyr::filter( is.na(param.smooth.df) | (param.smooth.df == "3L")) %>%
                dplyr::filter( !method == "unadjusted") %>%
                dplyr::filter( !(method %in% NAmethods)) %>%
                select( method, value ) %>%
                dplyr::rename( nrejects = value) %>%
                mutate( method = gsub("-df03", "", method)) %>%
                mutate( method = gsub("(-a)(.*)", "", method)) %>%
                na.omit() %>%
                mutate( rank = rank(nrejects) / n(),
                       propMaxRej = ifelse(nrejects == min(nrejects) & 
                                           nrejects == max(nrejects) &
                                           nrejects == 0, 0, 
                                           nrejects / max(nrejects, na.rm = TRUE)),
                       propPossible = nrejects / nrow(x)) %>%
                mutate(casestudy = colLabels[i])
            
            ranks <- rbind(ranks, tmp)
        }
    }
    return(ranks)
}

