### load packages
library(dplyr) # select part of data
library(stringr) # count separators
library(data.table) # transpose data frame
library(impute)
library(imputeLCMD)
library(tibble) # data frame manipulation

library(tidyverse)
library(ggplot2)#, # plots
library(viridis) # colorblind save color schemes
library(GGally) # paralell plot
library(fmsb) # spider chart
library(scales) # scale opacity of filling (alpha)
library(ggpubr) # multiple plots on one page
library(ggmosaic)


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

library(boot)

library(qgraph)


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
project <- "meat"

working_directory <- "/home/lisa/FH/Masterarbeit/LipidomeComparison"
setwd(working_directory)

lipid_list_path <- "/home/lisa/FH/Masterarbeit/LipidomeComparison/data/meat_fish_final_raw.csv"
annotation_path <- "/home/lisa/FH/Masterarbeit/LipidomeComparison/data/meat_annotation.csv"
data_matrix_path <- paste("/home/lisa/FH/Masterarbeit/LipidomeComparison/data/", project, "_data_matrix.csv", sep = "")

plot_path <- paste(working_directory, "/plots", sep = "")
plot_name <- paste(plot_path, "/lipid_list", sep = "")

## meat data
### data processing

lipid_list <- read.csv(lipid_list_path, sep = ",", dec = ".", header = TRUE)
annotation_data <- read.csv(annotation_path, sep = ",", dec = ".", header = TRUE)

lipid_list <- subset(lipid_list, select = c(Compound, Type, Filename, Status, Area))
lipid_list <- subset(lipid_list, Status == "Processed")
# lipid_list[lipid_list==''] <- NA
lipid_list[lipid_list=='N/F'] <- NA
lipid_list$Area <- as.numeric(as.character(lipid_list$Area))
target_lipids <- subset(lipid_list, Type == "Target Compound")
lipid_standards <- subset(lipid_list, Type == "Internal Standard")

target_lipids <- flip_df(target_lipids)
lipid_standards <- flip_df(lipid_standards)

target_lipids$SID <- sub(".*probe","sample", target_lipids$SID)
target_lipids$SID <- sub("\\.*pos2","2", target_lipids$SID)
target_lipids$SID <- sub("\\.*pos","1", target_lipids$SID)

target_lipids <- subset(target_lipids, str_detect(target_lipids$SID, "sample") == TRUE)

meat_AS <- target_lipids$SID[str_detect(target_lipids$SID, "AS") == TRUE]
target_lipids$SID[str_detect(target_lipids$SID, "AS") == TRUE] <- sub(".*sample","AS_sample", meat_AS)
meat_N <- target_lipids$SID[str_detect(target_lipids$SID, "AS") == FALSE]
target_lipids$SID[str_detect(target_lipids$SID, "AS") == FALSE] <- sub(".*sample","N_sample", meat_N)
target_lipids$SID <- str_remove(target_lipids$SID, "_AS")

meta_info <- read.table(text = target_lipids$SID, sep = "_")
colnames(meta_info) <- c("Treatment", "Sample_nr", "Biol_rep", "Tech_rep")
meta_info$Biol_rep <- paste(meta_info$Sample_nr, meta_info$Biol_rep, sep = "_")
meta_info$Tech_rep <- paste(meta_info$Biol_rep, meta_info$Tech_rep, sep = "_")

target_lipids <- cbind(target_lipids$SID, meta_info, target_lipids[, -1])
target_lipids <- droplevels(target_lipids)

map <- data.frame(Sample_nr=c("sample1","sample2","sample3", "sample4", "sample5", "sample6", "sample7"), 
                  Group_new=c("beef", "beef", "beef", "game", "game", "fish", "beef"))
target_lipids <- left_join(target_lipids, map, by="Sample_nr")
target_lipids$Group <- target_lipids$Group_new
target_lipids <- target_lipids[-ncol(target_lipids)]

# levels(target_lipids$Group)[levels(target_lipids$Group) == "fleisch"] <- "meat"
# levels(target_lipids$Group)[levels(target_lipids$Group) == "wild"] <- "game"
# levels(target_lipids$Group)[levels(target_lipids$Group) == "FISCH"] <- "fish"
# colnames(target_lipids) <- c("SID", colnames(target_lipids[-1]))
# meat_N <- subset(target_lipids, Treatment == "N")
# meat_AS <- subset(target_lipids, Treatment == "AS")

## Exploratory data analysis

### log2
log_meat <- log2(select_if(target_lipids, is.numeric))
log_meat <- cbind(select_if(target_lipids, is.character), select_if(target_lipids, is.factor), log_meat)
lipid_list <- log_meat

#### remove columns where all values are missing
impute_meat <- lipid_list[, which(colMeans(!is.na(lipid_list)) > 0.8)]
impute_meat <- as.matrix(select_if(impute_meat, is.numeric))

#### perform missing data imputation
meat_imputed <- as.data.frame(impute.QRILC(impute_meat, tune.sigma = 1)[[1]])
# meat_imputed <- as.data.frame(impute.MinDet(impute_meat))
meat_imputed <- cbind(lipid_list[, 1:6], meat_imputed)
meat_imputed <- droplevels(meat_imputed) # remove unused levels from factors

#### calculate the means for the replicates
meat_groups <- generate_categorical_table(meat_imputed$Group)

meat_numeric <- meat_imputed
# meat_numeric$Sample_nr <- as.numeric(meat_numeric$Sample_nr)
meat_numeric$Group <- as.numeric(meat_numeric$Group)
# meat_by_sample <- calc_by_replicate(meat_numeric, meat_numeric$Sample_nr, mean)
meat_by_replicate <- calc_by_replicate(meat_numeric, meat_numeric$Biol_rep, mean)

lipid_list <- paste_catecorical_variable(meat_by_replicate, ncol(meat_by_replicate), meat_groups)

                                         
colnames(lipid_list)[1] <- "Bio_replicate"
# lipid_list$Sample_nr <- paste("sample", lipid_list$Sample_nr, sep = "")

### graphical exploratory data analysis
# qqplot_by_factor(lipid_list, "Group", out_path = plot_name)
# histogram_by_factor(lipid_list, "Group", out_path = plot_name)
# boxplot_by_factor(lipid_list, "Group", out_path = plot_name)
#
# parallel_plot(lipid_list, lipid_list$Group, out_path = plot_name)
# meat_spider <- calc_by_replicate(lipid_list, lipid_list$Group, mean)
# spider_chart(meat_spider, legend_lab = meat_spider$Group.1, out_path = plot_name)


########### Lipid information #############

lipids <- colnames(select_if(lipid_list, is.numeric))
lipid_class <- sub("\\ .*","", lipids)

my_sub <- sub(".*\\(", "", lipids)
total_chain <- sub("\\-.*", "", my_sub)
total_chain <- sub("\\).*", "", total_chain)
fa <- sub(".*\\-", "", my_sub)
fa <- sub("\\).*", "", fa)
fa <- sub(".*e.*", NA, fa)
fa <- ifelse(is.na(fa), total_chain, fa)

lipid_species <- paste(head_group, total_chain)
fatty_acyls <- paste(head_group, fa)

saturated_bonds <- as.numeric(sub(".*\\:", "", lipid_meta$lipidSpecies))
is_pufa <- saturated_bonds >= 2

lipid_meta <- data.frame(lipidID = lipids, 
                         lipidClass = lipid_class, 
                         lipidSpecies = total_chain, 
                         fattyAcyls = fa, 
                         PUFA = is_pufa)

class_table <- as.data.frame(table(lipid_class))
species_table <- as.data.frame(table(sub(".*\\ ", "", lipid_meta$lipidSpecies)))
alcyl_table <- as.data.frame(table(sub(".*\\ ", "", lipid_meta$fattyAcyls)))

class_plot <- ggplot(data = lipid_meta) +
  geom_mosaic(aes(x = product(lipidClass), fill=lipidClass)) +
  scale_fill_viridis_d() +
  labs(x="Lipid class", y = NULL, title='Frequency of lipid classes')

pufa_plot <- ggplot(data = lipid_meta) +
  geom_mosaic(aes(x = product(PUFA), fill=PUFA)) +
  scale_fill_viridis_d(name = "", labels = c("no PUFA", "PUFA")) +
  labs(x="PUFA", y = NULL, title='Proportion of lipids containing at least one PUFA') + 
  theme(legend.position = "bottom", legend.title = element_blank())





# x <- layer_data(p, 1)
# y <- mutate(x , m.x = (xmin + xmax)/2, m.y =  (ymin + ymax)/2)
# z <- subset(y, select = c(m.x, m.y))
# 
# p + geom_text_repel(
#     # extract rectangle centers, add labels
#     data = z,
#     aes(x = m.x, y = m.y, label = levels(lipid_meta$lipidClass), 
#         angle = 90)
#   )
  
  

    
### test for normality
meat_normality <- shapiro_by_factor(lipid_list, lipid_list$Group)

### test for correlation
meat_correlation <- cor(select_if(lipid_list, is.numeric), method = "spearman")
matrix_heatmap(meat_correlation, title = "Correlation heatmap", interactive = FALSE)

## correlation network plot
library(corrr)
meat_cor <- cor(select_if(lipid_list, is.numeric))
meat_cor
meat_cor %>% fashion()
select_if(lipid_list, is.numeric) %>% correlate() %>%
  network_plot(min_cor = 1)



corMat <- cor(select_if(lipid_list, is.numeric)) # Correlate data
qgraph(corMat, graph = "pcor", layout = "spring", theme = "TeamFortress", 
       posCol = viridis(n = 1), negCol = viridis(n = 1, direction = -1))



### lipid ratios
meat_sample_means <- calc_by_replicate(lipid_list, factor = lipid_list$Sample_nr, mean)
meat_group_means <- calc_by_replicate(lipid_list, factor = lipid_list$Group, mean)

## ratio heatmap
{
# ratio matrix
ratio_list <- lapply(1:nrow(meat_group_means),
                       function(i) calculate_ratio_matrix(as.numeric(meat_group_means[i, -1]),
                                                          names_vector = colnames(meat_group_means[])[-1]))
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
beef <- subset(lipid_list, subset = lipid_list$Group == "beef")
beef$Sample_nr <- substr(beef$Bio_replicate, start = 1, stop = 7)
map <- data.frame(Sample_nr=as.factor(annotation_data$Sample), 
                  Nutrition=annotation_data$Nutrition)
beef <- left_join(beef, map, by="Sample_nr")
beef$Nutrition <- droplevels(beef$Nutrition)
beef$Sample_nr <- as.factor(beef$Sample_nr)


beef_pca <- PCA(select_if(beef, is.numeric), scale.unit = TRUE, graph = FALSE)
beef_eig <- as.data.frame(get_eigenvalue(beef_pca))
scree_factoextra(beef_pca)
scree_base(select_if(beef, is.numeric))
bp <- biplot_ggplot2(beef, "Nutrition", loadings = FALSE, ellipse = T, scale = T)
bp + geom_text(aes(label = substr(beef$Sample_nr, start = 7, stop = 7), hjust = 1))
biplot_factoextra(beef_pca, beef$Nutrition, ellipse = TRUE)

fviz_pca_ind(beef_pca, 
                habillage = beef$Nutrition, # a vector of groups by whicht to color   
                label = "var", 
                labelsize = 2, 
                geom.var = "none") + 
  geom_text(aes(breaks = beef$Sample_nr, label = substr(beef$Sample_nr, start = 7, stop = 7), hjust = 1))
  

meat_pca <- PCA(select_if(lipid_list, is.numeric), scale.unit = TRUE, graph = FALSE)
meat_eigenvalue <- as.data.frame(get_eigenvalue(meat_pca))
scree_factoextra(meat_pca)
simple_barchart(data_frame = meat_eigenvalue, y = meat_eigenvalue$eigenvalue, x = c(1:nrow(meat_eigenvalue)))
scree_base(select_if(lipid_list, is.numeric))
biplot_ggplot2(lipid_list, "Group", loadings = FALSE, ellipse = TRUE, scale = TRUE)
biplot_factoextra(meat_pca, lipid_list$Group, ellipse = TRUE)

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

# most contributional lipids
meat_var <- meat_pca$var
meat_contrib <- as.data.frame(meat_var$contrib)
pc1_contrib_table <- meat_contrib[order(meat_contrib$Dim.1, decreasing = TRUE),]
pc1_contrib <- rownames(pc1_contrib_table)[1:10]
pc2_contrib_table <- meat_contrib[order(meat_contrib$Dim.2, decreasing = TRUE),]
pc2_contrib <- rownames(pc2_contrib_table)[1:10]
meat_pca_sub <- subset(lipid_list, select = c("Sample_nr", "Group", pc1_contrib, pc2_contrib))

parallel_plot(meat_pca_sub, meat_pca_sub$Group)

means_meat_pca <- calc_by_replicate(meat_pca_sub, meat_pca_sub$Group, funct = mean)
# spider_chart(means_meat_pca, legend_lab = means_meat_pca$)

### PLS-DA



### Clustering
meat_clust <- data.frame(Group = lipid_list$Group)
meat_clust <- cbind(meat_clust, select_if(lipid_list, is.numeric))
rownames(meat_clust) <- lipid_list$SID

hclust_performance_table(meat_clust)
hclust_performance_plot(meat_clust)

t_lipid_list <- as.data.frame(t(select_if(lipid_list, is.numeric)))
colnames(t_lipid_list) <- lipid_list$Bio_replicate

hclust_performance_plot(t_lipid_list)
lipid_hclust <- hclust(dist(t_lipid_list, method = "manhattan"), method = "average")
hclust_dendrogram(lipid_hclust)

meat_dist <- dist(select_if(meat_clust, is.numeric), method = "manhattan")
meat_hclust <- hclust(meat_dist, method = "average")
hclust_dendrogram(meat_hclust,
                  labs = paste(lipid_list$Sample_nr,
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
meat_kruskal <- kruskal_test_by_col(lipid_list, "Group")
meat_kruskal$p_adj <- p.adjust(meat_kruskal$p_value, method = "fdr")
meat_significant_k <- subset(meat_kruskal, meat_kruskal$p_value <= 0.05)

# anova
meat_anova <- one_way_anova_by_col(lipid_list, "Group")
meat_anova$p_adj <- p.adjust(meat_anova$p_value, method = "fdr")
meat_significant_a <- subset(meat_anova, meat_anova$p_value <= 0.05)


{ # meat vs fish volcano plot
meat_vs_fish <- subset(lipid_list, Group == "fish" | Group == "beef")
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
             significance = 0.05)
}

{ # meat vs game volcano plot
  meat_vs_game <- subset(lipid_list, Group == "game" | Group == "beef")
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
               foldchange = 2)
}

{ # game vs fish volcano plot
  # game vs fish volcano plot
  game_vs_fish <- subset(lipid_list, Group == "fish" | Group == "game")
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

