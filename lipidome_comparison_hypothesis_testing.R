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

working_data

p_values_by_column <- function(input_df,
                           test_method = t.test,
                           alternative = "two.sided"){
  
  p_values <- apply(select_if(input_df, is.numeric), 2, 
                    function(x) test_method(x ~ input_df$treatment, alternative = alternative)$p.value)
  as.data.frame(p_values)
}

log2_foldchange <- function(input_df){
  log2_df <- log2(select_if(input_df, is.numeric))
  log2_df$treatment <- input_df$treatment
  
  means <- aggregate(select_if(log2_df, is.numeric), by = list(log2_df$treatment), FUN = mean)
  rownames(means) <- means$Group.1
  means <- as.data.frame(select_if(means, is.numeric))
  log2_foldchange <- vector()
  for(i in 1:ncol(means)){
    log2_foldchange[i] <- means[1, i] - means[2, i]
  }
  
  as.data.frame(log2_foldchange)
}

p_vals <- p_values_by_column(working_data)
adj_p_vals <- p.adjust(p_vals$p_values, method = "fdr")
log2FC <- log2_foldchange(working_data)

volcano_df <- cbind(p_vals, adj_p_vals, log2FC)

threshold_p <- volcano_df$p_values < 0.05
threshold_adj <- volcano_df$adj_p_vals < 0.05
new_volcano_df <- cbind(volcano_df, threshold_p, threshold_adj)


volcano_plot <- function(volcano_df, 
                         x, y, 
                         significance = 0.05,
                         foldchange = 0.05, 
                         color_by = "significance",
                         title = "Volcano plot", 
                         x_lab = "log2Fold", y_lab = "-log10(p-value)"){

  if(color_by == "significance"){
    threshold_col <- y < significance
    threshold_shape <- abs(x) < foldchange
  } 
  else if(color_by == "foldchange"){
    threshold_col <- abs(x) < foldchange
    threshold_shape <- y < significance
  }

  print(threshold_shape)
  limits <- max(-1*min(x), max(x))
  
  new_volcano_df <- cbind(volcano_df, threshold_col)
  
  volcano <- ggplot(data = new_volcano_df, 
                    aes(x = x, y = -1*log10(y), 
                        col = threshold_col
                        # shape = threshold_shape
                        )) + 
    geom_point() + 
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
    scale_color_viridis_d(begin = 0, 
                          end = 0.8, 
                          labels = c("not significant", "significant"), 
                          name = "Significance")

  volcano
}

volcano_plot(new_volcano_df, new_volcano_df$log2_foldchange, new_volcano_df$p_values, 0.05)



