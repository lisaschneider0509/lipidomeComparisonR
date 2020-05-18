### load packages
library(dplyr) # select part of data
library(stringr) # count separators
library(data.table) # transpose data frame
library(impute)
library(imputeLCMD)

library(ggplot2)#, # plots
library(tibble) # data frame manipulation
library(viridis) # colorblind save color schemes
library(GGally) # paralell plot
library(fmsb) # spider chart
library(scales) # scale opacity of filling (alpha)
library(ggpubr) # multiple plots on one page

library(ggrepel)
library(factoextra)
library(ggfortify) # biplot with ggplot
library(corrplot)
library(FactoMineR)

library(heatmaply) # interactive heatmap
library(gplots) # heatmap
library(plotly) # interactive ggplots
library(htmlwidgets) # save plotly-plots as html
library(dendextend)


# library(psych) # for correlation plot 
# library(gridExtra)
# library(devtools)
# library(dendextend)
# library(limma) # hypothesis testing


source("lipidome_comparison_dataTransformaions.R")
source("lipidome_comparison_EDA.R")
source("lipidome_comparison_pca.R")
source("lipidome_comparison_clustering.R")
source("lipidome_comparison_hypothesis_testing.R")

# set ggplot theme
my_theme <- theme_set(
  theme_minimal() +
    theme(plot.title = element_text(size=12, hjust = 0.5, family="AvantGarde"),
          plot.subtitle = element_text(size = 8, hjust = 0.5, family = "AvantGarde", colour = "grey40"),
          axis.text.x = element_text(size = 8, colour = "grey40", family="AvantGarde"),
          axis.text.y = element_text(size = 8, colour = "grey40", family="AvantGarde"),
          axis.title.x = element_text(size = 10, hjust = 0.5, colour = "grey40", family="AvantGarde"),
          axis.title.y = element_text(size = 10, hjust = 0.5, colour = "grey40", family="AvantGarde"),
          legend.text = element_text(size = 8, colour = "grey40", family="AvantGarde"),
          legend.title = element_text(size = 10, colour = "grey40", family="AvantGarde"))
)

## set variables
working_directory <- "/home/lisa/FH/Masterarbeit/LipidomeComparison"
setwd(working_directory)

test_path <- "/home/lisa/FH/Masterarbeit/LipidomeComparison/data/Probe-Datensatz_lisa.csv"
meat_data_path <- "/home/lisa/FH/Masterarbeit/LipidomeComparison/data/meat_fish_final_raw.csv"

plot_path <- paste(working_directory, "/plots", sep = "")
plot_name <- paste(plot_path, "/meat_data", sep = "")

## meat data
### data processing

meat_data <- read.csv(meat_data_path, sep = ",", dec = ".", header = TRUE)
meat_data <- subset(meat_data, select = c(Compound, Type, Filename, Status, Group, Area))
meat_data <- subset(meat_data, Status == "Processed")
# meat_data[meat_data==''] <- NA
meat_data[meat_data=='N/F'] <- NA
meat_data$Area <- as.numeric(as.character(meat_data$Area))
meat_target <- subset(meat_data, Type == "Target Compound")
meat_standard <- subset(meat_data, Type == "Internal Standard")

meat_target <- flip_df(meat_target)

meat_target$SID <- sub(".*probe","sample", meat_target$SID)
meat_target$SID <- sub("\\.*pos2","2", meat_target$SID)
meat_target$SID <- sub("\\.*pos","1", meat_target$SID)

meat_target <- subset(meat_target, str_detect(meat_target$SID, "sample") == TRUE)

meat_AS <- meat_target$SID[str_detect(meat_target$SID, "AS") == TRUE]
meat_target$SID[str_detect(meat_target$SID, "AS") == TRUE] <- sub(".*sample","AS_sample", meat_AS)
meat_N <- meat_target$SID[str_detect(meat_target$SID, "AS") == FALSE]
meat_target$SID[str_detect(meat_target$SID, "AS") == FALSE] <- sub(".*sample","N_sample", meat_N)
meat_target$SID <- str_remove(meat_target$SID, "_AS")

meta_info <- read.table(text = meat_target$SID, sep = "_")
colnames(meta_info) <- c("Treatment", "Sample_nr", "Biol_rep", "Tech_rep")
meta_info$Biol_rep <- paste(meta_info$Sample_nr, meta_info$Biol_rep, sep = "_")
meta_info$Tech_rep <- paste(meta_info$Biol_rep, meta_info$Tech_rep, sep = "_")

meat_target <- cbind(meat_target$SID, meta_info, meat_target[, -1])
meat_target <- droplevels(meat_target)

map <- data.frame(Sample_nr=c("sample1","sample2","sample3", "sample4", "sample5", "sample6", "sample7"), 
                  Group_new=c("beef", "beef", "beef", "game", "game", "fish", "beef"))
meat_target <- left_join(meat_target, map, by="Sample_nr")
meat_target$Group <- meat_target$Group_new
meat_target <- meat_target[-ncol(meat_target)]

# levels(meat_target$Group)[levels(meat_target$Group) == "fleisch"] <- "meat"
# levels(meat_target$Group)[levels(meat_target$Group) == "wild"] <- "game"
# levels(meat_target$Group)[levels(meat_target$Group) == "FISCH"] <- "fish"
# colnames(meat_target) <- c("SID", colnames(meat_target[-1]))
# meat_N <- subset(meat_target, Treatment == "N")
# meat_AS <- subset(meat_target, Treatment == "AS")

meat_data <- meat_target

## Exploratory data analysis

### impute missing values #todo find option to avoid imputation of negative values
#### remove columns where all values are missing
impute_meat <- meat_data[, which(colMeans(!is.na(meat_data)) > 0.8)] 
impute_meat <- as.matrix(select_if(impute_meat, is.numeric))

#### perform missing data imputation
# meat_imputed <- as.data.frame(impute.QRILC(impute_meat, tune.sigma = 1)[[1]])
meat_imputed <- as.data.frame(impute.MinDet(impute_meat))
meat_imputed <- cbind(meat_data[, 1:6], meat_imputed)
meat_imputed <- droplevels(meat_imputed) # remove unused levels from factors

#### calculate the means for the replicates
meat_groups <- generate_categorical_table(meat_imputed$Group)

meat_numeric <- meat_imputed
meat_numeric$Sample_nr <- as.numeric(meat_numeric$Sample_nr)
meat_numeric$Group <- as.numeric(meat_numeric$Group)
# meat_by_sample <- calc_by_replicate(meat_numeric, meat_numeric$Sample_nr, mean)
meat_by_replicate <- calc_by_replicate(meat_numeric, meat_numeric$Biol_rep, mean)

meat_data <- paste_catecorical_variable(meat_by_replicate, 3, meat_groups)
colnames(meat_data)[1] <- "Bio_replicate"
meat_data$Sample_nr <- paste("sample", meat_data$Sample_nr, sep = "")

### graphical exploratory data analysis
# qqplot_by_factor(meat_data, "Group", out_path = plot_name)
# histogram_by_factor(meat_data, "Group", out_path = plot_name)
# boxplot_by_factor(meat_data, "Group", out_path = plot_name)
# 
# parallel_plot(meat_data, meat_data$Group, out_path = plot_name)
# meat_spider <- calc_by_replicate(meat_data, meat_data$Group, mean)
# spider_chart(meat_spider, legend_lab = meat_spider$Group.1, out_path = plot_name)

### test for normality
meat_normality <- shapiro_by_factor(meat_data, meat_data$Group)

### test for correlation
meat_correlation <- cor(select_if(meat_data, is.numeric), method = "spearman")
matrix_heatmap(meat_correlation, title = "Correlation heatmap", interactive = FALSE)

### lipid ratios
meat_sample_means <- calc_by_replicate(meat_data, factor = meat_data$Sample_nr, mean)
meat_group_means <- calc_by_replicate(meat_data, factor = meat_data$Group, mean)

## ratio heatmap
{
# ratio matrix
ratio_list <- lapply(1:nrow(meat_group_means), 
                       function(i) calculate_ratio_matrix(as.numeric(meat_group_means[i, -1]), 
                                                          names_vector = colnames(meat_group_means[-1])))
names(ratio_list) <- as.character(meat_group_means$Group.1)


# save ratio matrixes to csv
lapply(1:length(ratio_list), 
       function(i) write.csv(ratio_list[[i]], 
                             file = paste("output/ratio_matrix_", 
                                          names(ratio_list[i]), ".csv", 
                                          sep = ""),
                             row.names = TRUE))
# ratio heatmap
heatmap_list <- lapply(1:length(ratio_list), 
       function(i) matrix_heatmap(ratio_list[[i]], 
                                  title = names(ratio_list)[i], 
                                  interactive = TRUE))
names(heatmap_list) <- names(ratio_list)

# ratio heatmap html
lapply(1:length(heatmap_list), 
       function(i) saveWidget(heatmap_list[[i]], 
                              file = paste(plot_path, 
                                           "/ratio_heatmap_", 
                                           names(heatmap_list)[i], 
                                           ".html", sep = "")))
# ratio heatmaps arranged
arranged_ratio_heatmaps <- ggarrange(plotlist = heatmap_list, 
          ncol = 2, 
          nrow = 4, 
          align = "h", 
          widths = c(0.9, 0.9),
          common.legend = TRUE, 
          legend = "right")
arranged_ratio_heatmaps <- annotate_figure(arranged_ratio_heatmaps, 
                                top = text_grob("Lipid ratios per sample", 
                                                size = 16, family = "AvantGarde"))

# save ratio heatmaps
ggsave(filename = paste(plot_name, "/ratio_heatmap.png", sep = ""), 
       plot = arranged_ratio_heatmaps, width = 8, height = 12)
}

# ratio barplots
ratio_lipids <- colnames(meat_group_means[c(1, 4:6)])
barchart_list <- plot_ratio_barcharts(meat_group_means, ratio_lipids, ratio_lipids[1])
arranged_barcharts <- ggarrange(plotlist = barchart_list,
                                align = "h", 
                                widths = c(0.9, 0.9),
                                common.legend = TRUE, 
                                legend = "none")


### PCA
meat_pca <- PCA(select_if(meat_data, is.numeric), scale.unit = TRUE, graph = FALSE)
meat_eigenvalue <- get_eigenvalue(meat_pca)
scree_factoextra(meat_pca)
scree_base(select_if(meat_data, is.numeric))
biplot_ggplot2(meat_data, "Group", loadings = FALSE, ellipse = TRUE, scale = TRUE)
biplot_factoextra(meat_pca, meat_data$Group, ellipse = TRUE)

# loadings plots 
fviz_pca_var(meat_pca, # factoextra
             geom = c("point"), 
             col.var = "contrib", 
             gradient.cols = viridis(n = 3, direction = -1), 
             repel = TRUE)

plot_loadings(meat_pca, colour = TRUE, top_loadings = 10) # diy wth ggplot

# contribution to PCs
plot_contrib_to_pc(meat_pca)

fviz_contrib(meat_pca, choice = "var", axes = 1, top = 10, 
             fill = viridis(n = 1, begin = 0.3), color = viridis(n = 1, begin = 0.3), 
             ggtheme = my_theme)
fviz_contrib(meat_pca, choice = "var", axes = 2, top = 10, 
             fill = viridis(n = 1, begin = 0.3), color = viridis(n = 1, begin = 0.3), 
             ggtheme = my_theme, 
             linecolor = "black")


### Clustering 
meat_clust <- data.frame(Group = meat_data$Group)
meat_clust <- cbind(meat_clust, select_if(meat_data, is.numeric))
rownames(meat_clust) <- meat_data$SID

hclust_performance_table(meat_clust)
hclust_performance_plot(meat_clust)

meat_dist <- dist(select_if(meat_clust, is.numeric), method = "manhattan")
meat_hclust <- hclust(meat_dist, method = "average")
hclust_dendrogram(meat_hclust, 
                  labs = paste(meat_data$Sample_nr, 
                               meat_clust$Group, sep = "-"))

hclust_heatmap(meat_clust,
               dist_method = "manhattan",
               hclust_method = "average",
               row_names = meat_clust$Group)
hclust_heatmap_interactive(meat_clust, 
                           dist_method = "manhattan", 
                           hclust_method = "average")

### hypothesis testing & volcano plot

# kruskal-wallis 
meat_kruskal <- kruskal_test_by_col(meat_data, "Group")
meat_kruskal$p_adj <- p.adjust(meat_kruskal$p_value, method = "fdr")
meat_significant_k <- subset(meat_kruskal, meat_kruskal$p_value <= 0.05)

# anova
meat_anova <- one_way_anova_by_col(meat_data, "Group")
meat_anova$p_adj <- p.adjust(meat_anova$p_value, method = "fdr")
meat_significant_a <- subset(meat_anova, meat_anova$p_value <= 0.05)


{ # meat vs fish volcano plot
meat_vs_fish <- subset(meat_data, Group == "fish" | Group == "beef")
meat_vs_fish <- droplevels(meat_vs_fish)

p_meat_vs_fish <- one_sample_test_by_col(meat_vs_fish, meat_vs_fish$Group, method = t.test)
adj_meat_vs_fish <- p.adjust(p_meat_vs_fish$p_values, method = "fdr")
fc_meat_vs_fish <- log2_foldchange(meat_vs_fish, 
                                   meat_vs_fish$Group, 
                                   control_group = "fish", 
                                   test_group = "beef")

meat_fish_volcano <- data.frame(p_value = p_meat_vs_fish, adj_p_value = adj_meat_vs_fish, log2_foldchange = fc_meat_vs_fish)
meat_fish_volcano <- meat_fish_volcano[complete.cases(meat_fish_volcano),]

volcano_plot(meat_fish_volcano, 
             foldchange_col = meat_fish_volcano$log2_foldchange, 
             significance_col = meat_fish_volcano$adj_p_value, 
             foldchange = 1, 
             significance = 0.05,
             out_path = plot_name)
}

{ # meat vs game volcano plot
  meat_vs_game <- subset(meat_data, Group == "game" | Group == "meat")
  meat_vs_game <- droplevels(meat_vs_game)
  
  p_meat_vs_game <- one_sample_test_by_col(meat_vs_game, meat_vs_game$Group, method = t.test)
  adj_meat_vs_game <- p.adjust(p_meat_vs_game$p_values, method = "fdr")
  fc_meat_vs_game <- log2_foldchange(meat_vs_game, 
                                     meat_vs_game$Group, 
                                     control_group = "meat", 
                                     test_group = "game")
  
  meat_game_volcano <- data.frame(p_value = p_meat_vs_game, adj_p_value = adj_meat_vs_game, log2_foldchange = fc_meat_vs_game)
  meat_game_volcano <- meat_game_volcano[complete.cases(meat_game_volcano),]
  
  volcano_plot(meat_game_volcano, 
               foldchange_col = meat_game_volcano$log2_foldchange, 
               significance_col = meat_game_volcano$adj_p_value, 
               foldchange = 1, 
               out_path = plot_name)
}

{ # game vs fish volcano plot
  # game vs fish volcano plot
  game_vs_fish <- subset(meat_data, Group == "fish" | Group == "game")
  game_vs_fish <- droplevels(game_vs_fish)
  
  p_game_vs_fish <- one_sample_test_by_col(game_vs_fish, game_vs_fish$Group, method = t.test)
  adj_game_vs_fish <- p.adjust(p_game_vs_fish$p_values, method = "fdr")
  fc_game_vs_fish <- log2_foldchange(game_vs_fish, 
                                     game_vs_fish$Group, 
                                     control_group = "fish", 
                                     test_group = "game")
  
  volcano_df <- data.frame(p_value = p_game_vs_fish, adj_p_value = adj_game_vs_fish, log2_foldchange = fc_game_vs_fish)
  volcano_df <- volcano_df[complete.cases(volcano_df),]
  
  volcano_plot(volcano_df, 
               foldchange_col = volcano_df$log2_foldchange, 
               significance_col = volcano_df$adj_p_value, 
               foldchange = 1, 
               significance = 0.05,
               title = "Fish vs. game",
               out_path = plot_name)
}
