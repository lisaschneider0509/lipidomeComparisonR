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
wd <- lipid_data
rownames(wd) <- wd$X
wd_con <- subset(working_data, working_data$treatment == "Con")
wd_lps <- subset(working_data, working_data$treatment == "LPS")

t.test(wd_con$`11-HETE_10.43`, wd_lps$`11-HETE_10.43`)

welch_ttest <- function(input_df, group1, group2){
  x = input_df[group1]
  y = input_df[group2]
  x = as.numeric(x)
  y = as.numeric(y)  
  results = t.test(x, y)
  results$p.value
}

p_values <- apply(select_if(wd, is.numeric), 1, welch_ttest, group1 = c(1:6), group2= c(7:12))
p_df <- as.data.frame(p_values)
p_fdr <- p.adjust(p_df$p_values, method = "fdr")

log2_wd <- log2(select_if(wd, is.numeric))
mean_lps <- apply(log2_wd[, 1:6], 1, mean)
mean_con <- apply(log2_wd[, 7:12], 1, mean)
log2_foldchange <- mean_con - mean_lps
log2_foldchange <- as.data.frame(log2_foldchange)

result_df <- cbind(log2_foldchange, p_df)
result_df$p_fdr <- p_fdr
rownames(result_df) <- rownames(wd)

threshold <- result_df$p_fdr < 0.05
result_df$threshold <- threshold


volcano <- ggplot(data = result_df, aes(x = log2_foldchange, y = -1*log10(p_fdr)), fill = color)
volcano + geom_point()

