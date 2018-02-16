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


