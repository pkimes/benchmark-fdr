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
#' @param alpha numeric value  or vector indicating which alpha value(s) to use for IHW 
#'  (need to be among of the values specified in the bench object). Default 0.10.
#' @param linePlot logical indicating whether to plot a line plot in lieu of 
#'  a heatmap  
#' @param excludeMethods character vector containing the names of methods to
#'  exclude from the heatmap. Default is to exclude fdrreg methods which are 
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
#' @param maxColor character indicating the color of the text annotation of the
#'  maximum value of the variable indicated in \code{annotate}. 
#' @author Keegan Korthauer  

plotMethodRanks <- function(objects, colLabels, alpha = 0.10, 
                            xlab = "Case Study",
                            fill = c("propMaxRejections", "meanRank", 
                                     "meanRejections",
                                     "FDR", "TPR", "TNR"),
                            linePlot = FALSE,
                            excludeMethods = c("fdrreg-t",
                                               "fdrreg-e"),
                            rowOrder = NULL,
                            tableOnly = FALSE,
                            Nlabel = TRUE,
                            annotate = "propPossible",
                            maxColor= "white"){
  fill <- match.arg(fill)
  if (is.list(readRDS(objects[1]))){
    ranks <- data.frame()
    
    for (l in seq_along(objects)){
      x <- readRDS(objects[l])
      tmp <- tidy_df(x, colLabels = rep(colLabels[l], length(x)), 
                     fill, annotate, alpha) %>%
        group_by(method, casestudy) %>%
        summarize(nrejects = mean(nrejects),
                  rank = mean(rank),
                  propMaxRej = mean(propMaxRej),
                  propPossible= mean(propPossible))
      ranks <- rbind(ranks, as.data.frame(tmp))
    }
    
  }else{
    ranks <- tidy_df(objects, colLabels, fill, annotate, alpha)
  }
  
  # exclude methods 
  if(!is.null(excludeMethods)){
    ranks <- ranks %>% dplyr::filter( !(method %in% excludeMethods))
  }
  
  # add NAs for missing combinations of methods and casestudies
  ranks <- ranks %>% complete(method, casestudy)
  
  # average over datasets within case study
  if (fill %in% c("FDR", "TPR", "TNR")){
    annot <- "nrejects"
  }else{
    annot <- annotate
  }
  
  if (fill == "propMaxRejections"){
    annot2 <- "mean_prop"
  }else if(fill == "meanRank"){
    annot2 <- "mean_rank"
  }else{
    annot2 <- "mean_nrej"
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
  
  # replace all but max value of annotate variable with NA (within casestudy)
  maxvals <- ranks_avg %>% 
    group_by(casestudy) %>%
    summarize(maxann = max(get(annot2), na.rm=TRUE))
  ranks_avg <- left_join(ranks_avg, maxvals, by = "casestudy")
  ranks_avg <- ranks_avg %>%
    mutate(topLayer = ifelse(get(annot2) == maxann, topLayer, NA))
  
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
        geom_tile() + 
        scale_fill_distiller("% Rejected\n(relative to max)",
                             palette = "Blues", direction = 1, limits = c(0, 1),
                             labels = scales::percent) +
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
        geom_tile() + 
        scale_fill_distiller("Mean Rank", palette = "Blues", direction = 1) +
        theme_bw() +
        xlab(xlab) +
        ylab("Method") 
      
      if(!Nlabel)
        Fig <- Fig +
          labs(fill = "Rank") 
    }else{
      # reorder rows
      if (is.null(rowOrder))
        ranks_avg$method <- factor(ranks_avg$method,
                                 levels=unique(ranks_avg$method[order(ranks_avg$meta_nrej)]))
      
      Fig <- ggplot(ranks_avg, aes(x = studyname, y = method, fill = mean_nrej)) + 
        geom_tile() + 
        scale_fill_distiller(paste0("Mean ", fill), palette = "Blues", direction = 1) +
        theme_bw() +
        xlab(xlab) +
        ylab("Method")
    }
    
    if (!is.null(annotate)){
      Fig <- Fig + 
        geom_text(data = ranks_avg, 
                  aes(label = ifelse(is.na(topLayer), "", 
                                     scales::percent(topLayer))),
                  color = maxColor)
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


#' Function to create a tidy data frame of performance metrics 
#' 
#' Takes as input a vector of filepaths, each one pointing to a SummarizedBenchmark
#' results object, labels for each object, as well as various annotation variables.
#' Used internally by plotMethodRanks, but can also be used on its own.
#' 
#' @param objects a character vector containing the full filepaths to 
#'  SummarizedBenchmark results objects.
#' @param colLabels a character vector containing the labels to use for the 
#'  heatmap columns (same length and order as `objects`). 
#' @param alpha numeric value  or vector indicating which alpha value(s) to use for IHW 
#'  (need to be among of the values specified in the bench object). Default 0.10.
#' @param fill character indicating the outcome variable of interest. Default
#'  is 'propMaxRejections'. For 'FDR' and 'TPR' sb objects need to contain 
#'  these performance metrics.
#' @param annotate character indicating what should be plotted as text labels
#'  (heatmap plot only). Default is proportion of total possible rejections 
#'  (propPossible). Could also be another variable (from \code{fill} choices) 
#'  or NULL (for no text label annotations).
#' @author Keegan Korthauer  
tidy_df <- function(objects, colLabels, fill, annotate, alpha){
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
            tmp <- data.frame()
            for (lev in alpha){
              df <- estimatePerformanceMetrics(x, alpha=lev)
              df$alpha <- lev
              tmp <- bind_rows(tmp, zeroRejectionsFDR(df, tidy=TRUE))  
            }
            
            if (fill %in% c("TPR", "FDR", "TNR") && !(fill %in% tmp$performanceMetric))
                stop(fill, " is not found in performanceMetrics")
            
            if (!is.null(annotate)){
                if (annotate %in% c("TPR", "FDR", "TNR") &&
                    !(annotate %in% tmp$performanceMetric))
                    stop(annotate, " is not found in performanceMetrics")
            }
            
            tmp <- tmp %>%
                dplyr::filter( performanceMetric == pmcol) %>%
                dplyr::rename( method = blabel) %>%
                dplyr::filter( is.na(param.alpha) | (param.alpha == alpha)) %>%
                dplyr::filter( is.na(param.smooth.df) | (param.smooth.df == "3L")) %>%
                dplyr::filter( !method == "unadjusted") %>%
                dplyr::filter( !(method %in% NAmethods)) %>%
                dplyr::select( method, value, alpha ) %>%
                dplyr::rename( nrejects = value) %>%
                mutate( method = gsub("-df03", "", method)) %>%
                mutate( method = gsub("(-a)(.*)", "", method)) %>%
                mutate( rank = rank(nrejects) / n(),
                       propMaxRej = ifelse(nrejects == min(nrejects) & 
                                           nrejects == max(nrejects) &
                                           nrejects == 0, 0, 
                                           nrejects / max(nrejects, na.rm = TRUE)),
                       propPossible = nrejects / nrow(x)) %>%
                mutate(casestudy = colLabels[i])
            if ( is.list(objects)) {
             tmp <- tmp %>%
               mutate(replicate = i)
            }
            
            ranks <- rbind(ranks, tmp)
        }
    }
    return(ranks)
}


# function to replace observations of NA FDR when nrejects = 0 to 0
# useful when averaging over simulation reps so we don't penalize a method
# that almost never rejects anything but we take the mean over all the 
# times it does
#' @param tidy_df a data frame output from estimatePerformanceMetrics with the 
#' tidy = FALSE option, where 
#' each row is an observation of a performance metric for a method. Assumes that
#' method is in the blabel column, the performance metric is in the performanceMetric
#' column, and that the value of the performance metric is in the value column.
#' @param tidy logical, whether to return the metrics table in tidy format or not
#' (analagous to the same argument in tidyUpMetrics function). Default is TRUE.
#' @author Keegan Korthauer
zeroRejectionsFDR <- function(df, tidy = TRUE){
  df$FDR <- ifelse(df$rejections == 0, 0, df$FDR)
  
  if(tidy){
    valueCols <- c("TPR", "FDR", "TNR", "FNR", "rejections")
    tidy_df <- gather(as.data.frame(df), keys = valueCols) %>%
      dplyr::rename(performanceMetric = key)
    return(tidy_df)
  }else{
    return(df)
  }
}

