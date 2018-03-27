#' Function to create boxplots of median effect sizes across simulation reps,
#' stratified by significance and DE status
#' 
#' Created for the yeast48-simulation analysis. Assumes that effect sizes are
#' in the rowData. For each method, calculates the median effect
#' size for genes detected as DE (at specified alpha), separately for true DE 
#' genes and null genes. The set of genes to be plotted is specified by the DE 
#' argument.
#' 
#' @param sbl a list of summarized benchmark objects, with each item representing
#' a simulation replicate
#' @param alpha the alpha cutoff
#' @param DE whether to plot the effect sizes for DE genes (if false, will plot
#' the effect sizes for non-DE genes)
#' @param covname character that specifies column name of rowData to use for 
#' y-axis (e.g. effect size, or ind_covariate)
#' 
#' @return a ggplot object
plotEffectSize <- function(sbl, alpha=0.05, DE=TRUE, covname){

summarize_one_item <- function(object, alpha){
  object <- object[,!( grepl("^ihw", as.character( colData( object )$blabel ) )
                       & colData( object )$param.alpha != alpha )]
  
  df <- as.data.frame(cbind(as.matrix(rowData(object)), 
                            1*(assays(object)[["qvalue"]]) < 0.05))
  colnames(df)[1] <- "truth" 
  df <- df %>%
    select(truth, covname, bonf, bh,
           qvalue, contains("ihw"), ashs, "bl-df03", lfdr, "scott-theoretical", 
           "scott-empirical") %>%
    rename("scott-theoretical"="scott-t") %>%
    rename("scott-empirical"="scott-e") %>%
    gather(method, significant, -truth, -covname) %>%
    group_by(method, truth, significant) %>%
    summarize(mean_effect_size=mean(abs(get(covname))),
              med_effect_size=median(abs(get(covname))),
              min_effect_size=min(abs(get(covname))),
              max_effect_size=max(abs(get(covname))),
              botQ=quantile(abs(get(covname)),0.25),
              topQ=quantile(abs(get(covname)), 0.75)) %>%
    na.omit()
  df$significant <- as.factor(df$significant)
  levels(df$significant)=c("Not Significant", "Significant")
  df$truth <- as.factor(df$truth)
  levels(df$truth)=c("Not DE", "DE")
  return(df)
}

df <- lapply(sbl, summarize_one_item, alpha=alpha)
df <- bind_rows(df, .id = "rep")
df <- as.tibble(df)

if (DE){
  df <- df %>% filter(truth=="DE")
}else{
  df <- df %>% filter(truth=="Not DE")
} 

p <- ggplot(df, aes(x = as.factor(method), y=med_effect_size)) +
  geom_boxplot() +
  facet_grid(significant~.) +
  xlab("Method") +
  ylab(paste0("Median ", covname))

if (DE){
  p <- p + ggtitle(paste0(covname, " of True DE genes by significance at alpha ", 
                       alpha))
}else{
  p <- p + ggtitle(paste0(covname, " of Non-DE genes by significance at alpha ", 
                       alpha))
} 

return(p)
}


#' Function to create lineplot of percentage significance for different covariate 
#' value bins.
#' 
#' Assumes that covariates are found 
#' in the rowData. For each method, divides the features into a specified number
#' of bins based on the value of the specified covariate, and 
#' calculates the proportion of features in that bin found
#' significant (at specified alpha). If the input object is a list of 
#' SummarizedBenchmark objects, this quantity is also averaged over the list
#' items (e.g. if each item is a simulation replicate). The function will
#' exclude any methods with missing values for rejections.
#' 
#' @param sbl a summarized benchmark object or a list of summarized benchmark 
#' objects (with each item representing
#' a simulation replicate)
#' @param alpha the alpha cutoff for determining significance.
#' @param covname character that specifies column name of rowData to use for 
#' y-axis (e.g. effect size, or ind_covariate)
#' @param nbins positive inter value indicating how many bins to divide the 
#' covariate values into for the x-axis (lineplot) or rows (heatmap). 
#' Defaults to 50.
#' @param trans character indicating a transformation to apply to the scale
#' of the y-axis (lineplot) or fill (heatmap). For example, "log1p" applies 
#' the log(x+1) transform. The default NULL means no transformation. 
#' @param linePlot logical whether to plot a line plot (smoothed density) 
#' (default is TRUE). If false, will plot a heatmap.
#' 
#' @return a ggplot object
covariateLinePlot <- function(sbl, alpha=0.05, nbins = 50, 
                             covname, trans = NULL, 
                             linePlot = TRUE){
  
  summarize_one_item <- function(object, alpha, nbins){
    object <- object[,!( grepl("^ihw", as.character( colData( object )$blabel ))
                         & colData( object )$param.alpha != alpha )]
    
    df <- as.data.frame(cbind(as.matrix(rowData(object)), 
                              1*(assays(object)[["qvalue"]]) < 0.05))
    colnames(df)[1] <- "truth" 
    df <- df %>%
      select(truth, covname, bonf, bh,
             qvalue, contains("ihw"), ashs, "bl-df03", lfdr, "scott-theoretical", 
             "scott-empirical") %>%
      dplyr::rename("scott-t" = "scott-theoretical") %>%
      dplyr::rename("scott-e" = "scott-empirical") %>%
      mutate(bin = ntile(abs(get(covname)), nbins)) %>%
      gather(method, significant, -covname, -bin) %>%
      group_by(method, bin) %>%
      summarize(nsig = sum(significant),
                tot = sum(!is.na(significant))) %>%
      dplyr::filter(method != 'truth') %>%
      na.omit()
    return(df)
  }
  
  if (is.list(sbl)){
    df <- lapply(sbl, summarize_one_item, alpha=alpha, nbins=nbins)
    df <- bind_rows(df, .id = "rep")
    df <- as.tibble(df) %>%
        mutate(prop = nsig / tot) %>%
        group_by(method, bin) %>%
        summarize(nsig = mean(prop)*100,
                  se = sd(prop * 100) / sqrt(n())) %>%
        na.omit()
  }else if ("SummarizedBenchmark" %in% class(sbl)){
    df <- summarize_one_item(sbl, alpha=alpha, nbins=nbins) %>%
      mutate(nsig = nsig/tot*100)
  }else{
    stop("Input object must be either a SummarizedBenchmark object, or ",
         "a list of SummarizedBenchmark objects.")
  }
  
  if (linePlot){
    p <- ggplot(df, aes(x = bin/nbins, y = nsig, color = method)) +
      geom_line(alpha = 3/4) +
      ylab("Mean % Significant") +
      scale_x_continuous(labels = scales::percent) +
      xlab(paste0(covname, " percentile")) +
      labs(color="Method") + 
      viridis::scale_color_viridis("Method", discrete = TRUE,
                                   guide = guide_legend(ncol = 2))
    
    if (!is.list(sbl)){
      p <- p + ylab("% Significant") 
    }else{
      p <- p + geom_errorbar(aes(ymin = nsig - se, ymax = nsig + se),
                    width=0.01, alpha=1/3)
    }
    
    if(!is.null(trans)){
      p <- p + scale_y_continuous(trans=trans)
    }
  }else{
    p <- ggplot(df, aes(x = as.factor(method), y = bin/nbins)) +
         geom_raster(aes(fill = nsig)) +
      xlab("Method") +
      scale_y_continuous(labels = scales::percent) +
      ylab(paste0(covname, " percentile")) +
      labs(fill="Mean % Significant")
    
    if (!is.list(sbl)){
      p <- p + labs(fill = "% Significant")
    }
    
    mpoint <- median(df$nsig)
    
    if(!is.null(trans)){
      p <- p + scale_fill_gradientn(colors=RColorBrewer::brewer.pal(8, "BuGn"), 
                                    trans=trans, na.value = "white")
    }else{
      p <- p + scale_fill_gradientn(colors=RColorBrewer::brewer.pal(8, "BuGn"))
    }
  }
  
  p <- p + ggtitle(paste0(covname, " by significance at alpha ", alpha))
  
  return(p)
}


