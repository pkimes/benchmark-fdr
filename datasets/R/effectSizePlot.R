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
    dplyr::select(truth, covname, bonf, bh,
                  qvalue, contains("ihw"), starts_with("ash"), "bl-df03", lfdr,
                  matches("fdrreg"))
  df <- df %>%
    tidyr::gather(method, significant, -truth, -covname) %>%
    dplyr::group_by(method, truth, significant) %>%
    dplyr::summarize(mean_effect_size=mean(abs(get(covname))),
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
#' @param metric string specifying which metric to plot as a function of the
#'        covariate specified by \code{covname}. Must be one of
#'        "discoveries" (default), "FDR", "TDR", or "TNR".
#' @param trans character indicating a transformation to apply to the scale
#' of the y-axis (lineplot) or fill (heatmap). For example, "log1p" applies 
#' the log(x+1) transform. The default NULL means no transformation. 
#' @param transx same as \code{trans}, but for x-axis. Only for lineplot.
#' @param linePlot logical whether to plot a line plot (smoothed density) 
#' (default is TRUE). If false, will plot a heatmap.
#' 
#' @return a ggplot object
covariateLinePlot <- function(sbl, alpha=0.05, nbins = 25, 
                              covname, metric = c("discoveries", "FDR", "TDR", "TNR") ,
                              trans = NULL, transx = NULL,
                              linePlot = TRUE, palette=candycols){
  metric <- match.arg(metric)
    
  summarize_one_item <- function(object, alpha, nbins){
    object <- object[,!( grepl("^ihw", as.character( colData( object )$blabel ))
                         & colData( object )$param.alpha != alpha )]
    
    df <- as.data.frame(cbind(as.matrix(rowData(object)), 
                              1*(assays(object)[["qvalue"]]) < 0.05))
    colnames(df)[1] <- "truth" 
    df <- df %>%
      dplyr::select(truth, covname, bonf, bh,
                    qvalue, contains("ihw"), contains("ash"),
                    "bl-df03", lfdr, 
                    contains("fdrreg"), contains("adapt")) %>%
      dplyr::mutate(bin = ntile(abs(get(covname)), nbins)) %>%
      tidyr::gather(method, significant, -covname, -bin, -truth) %>%
      dplyr::group_by(method, bin) %>%
      dplyr::summarize(nfp = sum(significant * (1-truth), na.rm = TRUE),
                       ntp = sum(significant * truth, na.rm = TRUE),
                       ntn = sum((1-significant) * (1-truth), na.rm = TRUE),
                       nsig = sum(significant, na.rm = TRUE),
                       tot = sum(!is.na(significant)))
  }
  
  if (is.list(sbl)){
    df <- lapply(sbl, summarize_one_item, alpha=alpha, nbins=nbins)
    df <- bind_rows(df, .id = "rep")
    df <- as_tibble(df)
    if (metric == "discoveries") {
        df <- dplyr::mutate(df, prop = nsig / tot)
    } else if (metric == "FDR") {
        df <- dplyr::mutate(df, prop = nfp / nsig)
    } else if (metric == "TDR") {
        df <- dplyr::mutate(df, prop = ntp / nsig)
    } else if (metric == "TNR") {
        df <- dplyr::mutate(df, prop = ntn / tot)
    } else {
        stop("metric was not a valid input.")
    }
    df <- dplyr::group_by(df, method, bin) %>%
        dplyr::summarize(nsig = mean(prop)*100,
                         se = sd(prop * 100) / sqrt(n())) 
  }else if ("SummarizedBenchmark" %in% class(sbl)){
    df <- summarize_one_item(sbl, alpha=alpha, nbins=nbins) %>%
      dplyr::mutate(nsig = nsig/tot*100)
  }else{
    stop("Input object must be either a SummarizedBenchmark object, or ",
         "a list of SummarizedBenchmark objects.")
  }
  
  # standardize method names and add color palette
  df <- df %>%
    dplyr::mutate(Method = gsub("-df03", "", method)) %>%
    dplyr::mutate(Method = gsub("(-a)(.*)", "", Method)) 
  df <- dplyr::left_join(df, palette, by="Method") 
  
  col <- as.character(df$col)
  names(col) <- as.character(df$Method)
  
  lty <- as.character(df$lty)
  names(lty) <- as.character(df$Method)

  if (metric == "discoveries") {
      ystr <- "Percent rejected"
  } else {
      ystr <- metric
  }

  if (linePlot){
    p <- ggplot(df, aes(x = (2*bin-1)/(2*nbins), y = nsig, color = Method)) +
      geom_line(alpha = 0.85, aes(linetype=Method)) +
      ylab(paste0("Mean ", ystr)) +
      scale_x_continuous(labels = scales::percent) +
      xlab(paste0(covname, " percentile")) +
      scale_color_manual(values = col) +
      scale_linetype_manual(values = lty) + 
      theme_classic() + 
      theme(axis.title = element_text(face="bold"),
            plot.title = element_text(face="bold"))
    
    if (!is.list(sbl)){
      p <- p + ylab(ystr)
    }else{
      p <- p + geom_errorbar(aes(ymin = nsig - se, ymax = nsig + se),
                    width=0.02, alpha=0.5)
    }
    
    if(!is.null(trans)){
      p <- p + scale_y_continuous(trans = trans, labels = function(x) paste0(x, "%"))
    }
    
    if(!is.null(transx)){
      p <- p + scale_x_continuous(trans=transx)
    }
  }else{
    p <- ggplot(df, aes(x = as.factor(Method), y = (2*bin-1)/(2*nbins))) +
         geom_raster(aes(fill = nsig)) +
      xlab("Method") +
      scale_y_continuous(labels = scales::percent) +
      ylab(paste0(covname, " percentile")) +
      labs(fill=paste0("Mean ", ystr))
    
    if (!is.list(sbl)){
      p <- p + labs(fill = ystr)
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


