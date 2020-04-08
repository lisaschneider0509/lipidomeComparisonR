## set ggplot theme
my_theme <- theme_set(
  theme_minimal() +
    theme(plot.title = element_text(size=12, hjust = 0.5, family="AvantGarde"),
          axis.text.x = element_text(size = 8, colour = "grey40", family="AvantGarde"),
          axis.text.y = element_text(size = 8, colour = "grey40", family="AvantGarde"),
          # axis.title = element_text(size = 10, colour = "grey40", family="AvantGarde"),
          axis.title.x = element_text(size = 10, hjust = 0.5, colour = "grey40", family="AvantGarde"),
          axis.title.y = element_text(size = 10, hjust = 0.5, colour = "grey40", family="AvantGarde"),
          legend.text = element_text(size = 8, colour = "grey40", family="AvantGarde"),
          legend.title = element_text(size = 10, colour = "grey40", family="AvantGarde"))
)

#' Calculate p-values by column
#' 
#' @description `p_values_by_column` takes a data frame and a grouping vector and returns a data frame with the p-values of a given method. 
#' @details This function takes a data frame ans a grouping vector and performs a given hypothesis test on all numeric variables. The results are returned in a new data frame. 
#' @param input_df data frame. Contains at least one numeric variable. 
#' @param group_vector vector. Vector containing the groups to be compared. group_vector = input_df$group (or other vector of length nrow(input_df))
#' @param alternative string. Specifies the alternative hypothesis. c("two.sided" (default), "greater", "less")
#' @param method function. Options: t.test (default), wilcox.test
#' @param confidence_level numeric. Confidence level of the interval. Default = 0.95
#' @example 
#' my_iris <- subset(x = iris, Species == "setosa" | Species == "versicolor")
#' one_sample_test_by_col(my_iris, my_iris$Species)
#' one_sample_test_by_col(my_iris, my_iris$Species, method = wilcox.test)
one_sample_test_by_col <- function(input_df,
                            group_vector,
                            alternative = "two.sided", 
                            method = t.test,
                            confidence_level = 0.95){
  
  p_values <- apply(select_if(input_df, is.numeric), 2, 
                    function(x) method(x ~ group_vector, 
                                       alternative = alternative, 
                                       conf.level = confidence_level)$p.value)
  as.data.frame(p_values)
}

#' Calculate log2 foldchange
#' 
#' @description `log2_foldchange` calculates the log2 foldchange for two groups of a given data frame. 
#' @details This function takes a data frame and a grouping vector and calculates the log2 foldchange between two groups. 
#' The results are returned in a new data frame. 
#' @param input_df data frame. Contains at least one numeric variable. 
#' @param group_vector vector. Vctor containing the grouping variabe. group_vector = input_df$group (or any vector of length nrow(input_df))
#' @example 
#' my_iris <- subset(x = iris, Species == "setosa" | Species == "versicolor")
#' log2_foldchange(my_iris, my_iris$Species)
log2_foldchange <- function(input_df, 
                            group_vector){
  log2_df <- log2(select_if(input_df, is.numeric))
  log2_df$group <- group_vector
  
  means <- aggregate(select_if(log2_df, is.numeric), by = list(log2_df$group), FUN = mean)
  rownames(means) <- means$Group.1
  means <- as.data.frame(select_if(means, is.numeric))
  log2_foldchange <- vector()
  for(i in 1:ncol(means)){
    log2_foldchange[i] <- means[1, i] - means[2, i]
  }
  
  out_df <- as.data.frame(log2_foldchange)
  rownames(out_df) <- colnames(select_if(input_df, is.numeric))
  out_df
}

#' Volcano plot
#' 
#' @description `volcano_plot` prints a volcano plot using a data frame containing p-values and foldchanges
#' @details This function takes a data frame with columns for p-values and (log2)foldchanges and returns a dotplot resembelling a volcano plot. 
#' It uses the given significance and foldchange threshods to mark significantly up- and down-regulated elements by color and by displaying
#' threshold lines. 
#' @param input_df data frame. 
#' @param foldchange_col vector. Column of input_df containing the foldchange values. Format: input_df$foldchange_col
#' @param significance_col vector. Column of input_df containing the p-values. Format: input_df$significance_col
#' @param significance numeric. Value for the significance threshold. Default = 0.05.
#' @param foldchange numeric. Vaue for the foldchange threshold. Default = #todo
#' @param title string. Main title. Default = "Volcano plot"
#' @param x_lab string. x-axix title. Default = "log2Fold"
#' @param y_lab string. y-axis title- Default = "-log10(p-value)"
#' @param labels vector of length nrow(volcano_df). Default = "none". If lables is set, points above the threshold are labelled. 
#' @param out_path string. Path to save volcano plot do png. Default = "none". 
#' If "none", the plot is either printed or saved to a variabe. 
#' @example 
#' set.seed(100)
#' fold_changes <- c(rnorm(2000, 0, 2))
#' pvalues <- runif(n=2000, min=1e-50, max=.1)
#' volcano_test <- as.data.frame(cbind(fold_changes, pvalues))
#' volcano_plot(volcano_test, 
#'              foldchange_col = volcano_test$fold_changes, 
#'              significance_col = volcano_test$pvalues, foldchange = 1)
#' volcano_plot(volcano_test, 
#'              foldchange_col = volcano_test$fold_changes, 
#'              significance_col = volcano_test$pvalues, 
#'              foldchange = 1, 
#'              labels =c(1:nrow(volcano_test)))
#' \dontrun
#' {dir.create(paste(getwd(), "/examples", sep = ""), showWarnings = FALSE)
#' dir <- paste(getwd(), "/examples/vp", sep = "")
#' volcano_plot(volcano_test, 
#'              foldchange_col = volcano_test$fold_changes, 
#'              significance_col = volcano_test$pvalues, 
#'              foldchange = 1, significance = 0.01, 
#'              out_path = dir)}
volcano_plot <- function(volcano_df,
                         foldchange_col, significance_col,
                         significance = 0.05,
                         foldchange = 2,
                         title = "Volcano plot",
                         x_lab = "log2Fold", y_lab = "-log10(p-value)", 
                         labels = vector(),
                         out_path = "none"){

  is_significant <- significance_col < significance
  
  threshold <- vector()
  for(i in 1:length(is_significant)){
    if(is_significant[i] == TRUE && foldchange_col[i] < -1*foldchange){
      buffer <- "down"
    }
    else if(is_significant[i] == TRUE && foldchange_col[i] > foldchange){
      buffer <- "up"
    }
    else{
      buffer <- "not_sig"
    }
    threshold <- c(threshold, buffer)
  }

  volcano_df <- cbind(volcano_df, threshold)

  limits <- max(-1*min(foldchange_col), max(foldchange_col))
  mycolors <- viridis(n = 2, begin = 0, end = 0.9)

  volcano <- ggplot(data = volcano_df, 
                    aes(x = foldchange_col, y = -1*log10(significance_col))) + 
    geom_point(aes(color = as.factor(threshold)), shape = 20) + 
    geom_hline(yintercept = -1*log10(significance), 
               linetype = "dashed", 
               colour = "grey40") +
    geom_vline(xintercept = -1*foldchange, 
               linetype = "dashed", 
               colour = "grey40") +
    geom_vline(xintercept = foldchange, 
               linetype = "dashed", 
               colour = "grey40") +
    labs(title = title) + 
    xlab(x_lab) + ylab(y_lab) + 
    scale_x_continuous(limits = c(-1*limits, limits)) +
    scale_color_manual(name = "Threshold",
                       values = c("up" = mycolors[1], "down" = mycolors[2], "not_sig" = "grey"), 
                       labels = c("Down-regulated", "Not significant FC", "Up-regulated")) + 
    theme(legend.position = "right")
  
  if (length(labels) != 0){
    volcano <- volcano + geom_text_repel(aes(label=ifelse(threshold != "not_sig", as.character(labels),'')),
                                   hjust=0, vjust=0,
                                   size = 3)
  }

  if(out_path != "none"){
    print(paste("Saving plot to ", out_path, "_volcano.png", sep = ""))
    ggsave(paste(out_path, "_volcano.png", sep = ""),
           plot = volcano)
  }
  volcano
}



