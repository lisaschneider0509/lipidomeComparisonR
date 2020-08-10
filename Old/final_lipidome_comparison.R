### load packages
{library(dplyr) # select part of data
library(stringr) # count separators
library(data.table) # transpose data frame
library(crmn)
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

source("lipidome_comparison_dataTransformaions.R")
source("lipidome_comparison_EDA.R")
source("lipidome_comparison_pca.R")
source("lipidome_comparison_clustering.R")
source("lipidome_comparison_hypothesis_testing.R")}

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

############### set variables #############################
{
  project <- "meat"

  working_directory <- "/home/lisa/FH/Masterarbeit/LipidomeComparison"
  setwd(working_directory)
  
  data_dir <- "/home/lisa/FH/Masterarbeit/LipidomeComparison/data"
  lipid_list_path <- "/home/lisa/FH/Masterarbeit/LipidomeComparison/data/meat_fish_final_raw.csv"
  annotation_path <- "/home/lisa/FH/Masterarbeit/LipidomeComparison/data/meat_annotation.csv"
  data_matrix_path <- paste("/home/lisa/FH/Masterarbeit/LipidomeComparison/data/", project, "_data_matrix.csv", sep = "")
  
  plot_path <- paste(working_directory, "/plots", sep = "")
}

############## improt TF data #############################
{
  lipid_list <- read.csv(lipid_list_path, sep = ",", dec = ".", header = TRUE)
  annotation_data <- read.csv(annotation_path, sep = ",", dec = ".", header = TRUE)
  
  lipid_list <- subset(lipid_list, select = c(Compound, RT, Type, Filename, Status, Area))
  lipid_list <- subset(lipid_list, Status == "Processed")
  # lipid_list[lipid_list==''] <- NA
  lipid_list[lipid_list=='N/F'] <- NA
  lipid_list$Area <- as.numeric(as.character(lipid_list$Area))
  
  retention_time <- subset(lipid_list, select = c(Compound, RT, Type))
  retention_time <- get_retention_time(retention_time, Compound = "Compound", RT = "RT")
  
  # target_lipids <- lipid_list
  target_lipids <- subset(lipid_list, Type == "Target Compound")
  lipid_standards <- subset(lipid_list, Type == "Internal Standard")
  
  target_lipids <- flip_df(target_lipids)
  lipid_standards <- flip_df(lipid_standards)
  
}

######## get sample information ############################
{
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

  map <- data.frame(Sample_nr=annotation_data$Sample,
                    Group_new=annotation_data$Meat.type)
  target_lipids <- left_join(target_lipids, map, by="Sample_nr")
  target_lipids$Group <- target_lipids$Group_new
  target_lipids <- target_lipids[-ncol(target_lipids)]

}

######### Data transformation #####################
{
  ### log2
  log_meat <- log2(select_if(target_lipids, is.numeric))
  log_meat <- cbind(select_if(target_lipids, is.character), select_if(target_lipids, is.factor), log_meat)
  lipid_list <- log_meat
  
  #### remove columns where all values are missing
  impute_meat <- lipid_list[, which(colMeans(!is.na(lipid_list)) > 0.8)]
  impute_meat <- as.matrix(select_if(impute_meat, is.numeric))
  norm_lipids <- impute_meat
  
  

  #### perform missing data imputation
  meat_imputed <- as.data.frame(impute.QRILC(impute_meat, tune.sigma = 1)[[1]])
  # meat_imputed <- as.data.frame(impute.MinDet(impute_meat))
  
  ### normalize
  # meat_imputed <- norm_lipids <- med_normalize(meat_imputed)
  
  meat_imputed <- cbind(lipid_list[, 1:6], meat_imputed)
  meat_imputed <- droplevels(meat_imputed) # remove unused levels from factors
  
  
  #### calculate the means for the replicates
  meat_groups <- generate_categorical_table(meat_imputed$Group)
  
  meat_numeric <- meat_imputed
  meat_numeric$Group <- as.numeric(meat_numeric$Group)
  meat_by_replicate <- calc_by_replicate(meat_numeric, meat_numeric$Biol_rep, mean)
  lipid_list <- paste_catecorical_variable(meat_by_replicate, ncol(meat_by_replicate), meat_groups)
  colnames(lipid_list)[1] <- "Bio_replicate"
  
  lipid_list$Sample <- sub("\\_.*", "", lipid_list$Bio_replicate)
}

write.csv(file = paste(data_dir, "/", project, "_data_imputed.csv", sep = ""), x = lipid_list)

########### Get lipid information #############
{
  lipids <- colnames(select_if(lipid_list, is.numeric))
  lipid_class <- sub("\\ .*","", lipids)
  lipid_class <- sub(".*\\-", "", lipid_class)
  
  my_sub <- sub(".*\\(", "", lipids)
  lipid_Species <- sub("\\-.*", "", my_sub)
  lipid_Species <- sub("\\).*", "", lipid_Species)
  fatty_acyls <- sub(".*\\-", "", my_sub)
  fatty_acyls <- sub("\\).*", "", fatty_acyls)
  fatty_acyls <- sub(".*e.*", NA, fatty_acyls)
  fatty_acyls <- ifelse(is.na(fatty_acyls), lipid_Species, fatty_acyls)
  alkenyl_ether <- grepl("P", fatty_acyls)
  
  lipid_species <- paste(lipid_class, lipid_Species)
  fatty_acyls <- paste(lipid_class, fatty_acyls)
  
  saturated_bonds <- as.numeric(sub(".*\\:", "", lipid_species))
  is_pufa <- saturated_bonds >= 2
  
  lipid_meta <- data.frame(lipidID = lipids, 
                           lipidClass = lipid_class, 
                           lipidSpecies = lipid_Species, 
                           fattyAcyls = fatty_acyls, 
                           PUFA = is_pufa, 
                           ether = alkenyl_ether)
  
  class_table <- as.data.frame(table(lipid_class))
  species_table <- as.data.frame(table(sub(".*\\ ", "", lipid_meta$lipidSpecies)))
  alcyl_table <- as.data.frame(table(sub(".*\\ ", "", lipid_meta$fattyAcyls)))
  table(is_pufa)
  
  class_plot <- ggplot(data = lipid_meta) +
    geom_mosaic(aes(x = product(lipidClass), fill=lipidClass)) +
    scale_fill_viridis_d(alpha = 1) +
    labs(x="Lipid class", y = NULL, title='Lipid classes') + 
    theme(axis.text = element_text(family = "AvantGarde"), 
          legend.position = "bottom", legend.title = element_blank())
  
  pufa_plot <- ggplot(data = lipid_meta) +
    geom_mosaic(aes(x = product(PUFA), fill=PUFA)) +
    scale_fill_viridis_d(name = "", labels = c("no PUFA", "PUFA")) +
    labs(x="PUFA", y = NULL, title='Lipids with at least one PUFA') + 
    theme(legend.position = "bottom", legend.title = element_blank())
  
  
  
  freq_plot <- ggarrange(plotlist = list(class_plot, pufa_plot),
            ncol = 2,
            nrow = 1,
            align = "h",
            widths = c(0.9, 0.9),
            common.legend = FALSE,
            legend = "bottom", 
            labels = c("A", "B"), 
            font.label = list(size = 10, 
                              color = "grey40", 
                              face = "plain", 
                              family = "AvantGarde"))
  
  ggsave(filename = paste(plot_path, "Lipid_frequencies.png", sep = "/"), 
         plot = freq_plot, 
         device = "png", 
         width = 10, 
         height = 5
         )
 
}

short_lipid_names <- paste(lipid_meta$fattyAcyls, " (", rowid(lipid_meta$fattyAcyls), ")", sep = "")
colnames(lipid_list) <- c("replicate", short_lipid_names, "group", "sample")

########################################################################

############ Univariate testing ########################################

## ANOVA
meat_anova <-list_anova_by_col(lipid_list, "group")
meat_anova <- meat_anova[-1]
meat_pval <- one_way_anova_by_col(lipid_list, "group")
meat_pval$p_adj <- p.adjust(meat_pval$p_value, method = "fdr")
meat_significant_a <- subset(meat_pval, meat_pval$p_adj <= 0.05)
lapply(meat_anova, TukeyHSD)

pca_anova <- list_anova_by_col(meat_pca_sub, "group")
pca_anova <- pca_anova[-1]
meat_p <- one_way_anova_by_col(meat_pca_sub, "group")
meat_adj <- p.adjust(meat_p$p_value, method = "fdr")
lapply(pca_anova, TukeyHSD)

########################################################################

###################### PCA #############################################
meat_pca <- PCA(select_if(lipid_list, is.numeric), scale.unit = T, graph = F)
meat_eigenvalue <- as.data.frame(get_eigenvalue(meat_pca))
scree <- scree_factoextra(meat_pca)

{ # scree
eigen_plot <- plot_pc_variance(meat_eigenvalue[1:10,], 
                               x = seq(1:10), 
                               y = meat_eigenvalue$eigenvalue[1:10], 
                               title = "Eigenvalue", 
                               ylab = "eigenvalues")
var_plot <- plot_pc_variance(meat_eigenvalue[1:10,], 
                            x = seq(1:10), 
                            y = meat_eigenvalue$variance.percent[1:10], 
                            title = "Variance [%]")
cum_var_plot <- plot_pc_variance(meat_eigenvalue[1:10,], 
                                 x = seq(1:10), 
                                 y = meat_eigenvalue$cumulative.variance.percent[1:10], 
                                 title = "Cummulative variance [%]", 
                                 ylab = NULL,
                                 hjust = 1)

scree <- ggarrange(plotlist = list(var_plot, cum_var_plot), 
          nrow = 1, 
          ncol = 2, 
          widths = c(1, 1), 
          labels = c("(A)", "(B)"), 
          font.label = list(size = 10, 
                            color = "grey40", 
                            face = "plain", 
                            family = "AvantGarde")
            )
ggsave(paste(plot_path, "scree_all.png", sep = "/"), scree, device = "png", 
       height = 5, width = 10)
}

meat_biplot <- biplot_ggplot2(lipid_list, groups = "group", loadings = FALSE, ellipse = TRUE, scale = TRUE, title = "Scores plot of meat data")
ggsave(paste(plot_path, "biplot_all.png", sep = "/"), meat_biplot, device = "png", 
       height = 5, width = 7)


# loadings plots
fviz_pca_var(meat_pca, # factoextra
             geom = c("point"),
             col.var = "contrib",
             gradient.cols = viridis(n = 3, direction = -1),
             repel = TRUE)


short_lipid_names <- paste(lipid_meta$fattyAcyls, " (", rowid(lipid_meta$fattyAcyls), ")", sep = "")
new_meat_pca <- meat_pca
rownames(new_meat_pca$var$coord) <-short_lipid_names

plot_loadings(new_meat_pca, colour = TRUE, top_loadings = 10, xlab = "PC1", ylab = "PC2")
ggsave(paste(plot_path, "loadings1_all.png", sep = "/"), 
       device = "png", 
       width = 7, 
       height =5)
plot_contrib_to_pc(new_meat_pca)
cpc1 <- fviz_contrib(new_meat_pca, choice = "var", axes = 1, top = 10,
             fill = viridis(n = 1, begin = 0.3), color = viridis(n = 1, begin = 0.3),
             title = "Contribution to PC1",
             ggtheme = my_theme)
cpc2 <- fviz_contrib(new_meat_pca, choice = "var", axes = 2, top = 10,
             fill = viridis(n = 1, begin = 0.3), color = viridis(n = 1, begin = 0.3),
             title = "Contribution to PC2",
             ggtheme = my_theme,
             linecolor = "black", 
             xtickslab.rt = 45)
contrib_pc <- ggarrange(plotlist = list(cpc1, cpc2), 
                        nrow = 1, 
                        ncol = 2, 
                        widths = c(1, 1), 
                        labels = c("A", "B"), 
                        font.label = list(size = 10, 
                                          color = "grey40", 
                                          face = "plain", 
                                          family = "AvantGarde"))
ggsave(paste(plot_path, "contrib_pc_all.png", sep = "/"), 
       contrib_pc, 
       device = "png", 
       width = 10, 
       height = 5)

# most contributional lipids
meat_var <- meat_pca$var
meat_contrib <- as.data.frame(meat_var$contrib)
pc1_contrib_table <- meat_contrib[order(meat_contrib$Dim.1, decreasing = TRUE),]
pc1_contrib <- rownames(pc1_contrib_table)[1:10]
pc2_contrib_table <- meat_contrib[order(meat_contrib$Dim.2, decreasing = TRUE),]
pc2_contrib <- rownames(pc2_contrib_table)[1:10]
meat_pca_sub <- subset(lipid_list, select = c("replicate", "group", pc1_contrib, pc2_contrib))

pp <- parallel_plot(meat_pca_sub, meat_pca_sub$group)
ggsave(filename = paste(plot_path, "paralellPlot.png", sep = "/"), plot = pp, 
       width = 7, 
       height = 5)

means_meat_pca <- calc_by_replicate(meat_pca_sub, meat_pca_sub$Group, funct = mean)

########## Beef PCA #################
beef <- subset(lipid_list, lipid_list$group == "beef")
map <- data_frame(sample = annotation_data$Sample, 
                  nutrition = annotation_data$Nutrition)
beef <- left_join(beef, map, by = "sample")
beef <- droplevels(beef)
beef_pca <- PCA(select_if(beef, is.numeric), scale.unit = T, graph = FALSE)
beef_eig <- as.data.frame(get_eigenvalue(beef_pca))
scree <- scree_factoextra(beef_pca)

{ # scree
  eigen_plot <- plot_pc_variance(beef_eig[1:10,], 
                                 x = seq(1:10), 
                                 y = beef_eig$eigenvalue[1:10], 
                                 title = "Eigenvalue", 
                                 ylab = "eigenvalues")
  var_plot <- plot_pc_variance(beef_eig[1:10,], 
                               x = seq(1:10), 
                               y = beef_eig$variance.percent[1:10], 
                               title = "Variance [%]")
  cum_var_plot <- plot_pc_variance(beef_eig[1:10,], 
                                   x = seq(1:10), 
                                   y = beef_eig$cumulative.variance.percent[1:10], 
                                   title = "Cummulative variance [%]", 
                                   ylab = NULL,
                                   hjust = 1)
  
  scree <- ggarrange(plotlist = list(var_plot, cum_var_plot), 
                     nrow = 1, 
                     ncol = 2, 
                     widths = c(1, 1), 
                     labels = c("(A)", "(B)"), 
                     font.label = list(size = 10, 
                                       color = "grey40", 
                                       face = "plain", 
                                       family = "AvantGarde")
  )
  ggsave(paste(plot_path, "scree_beef.png", sep = "/"), scree, device = "png", 
         height = 4, width = 10)
}

beef_biplot <- biplot_ggplot2(beef, groups = "nutrition", loadings = FALSE, ellipse = TRUE, scale = TRUE, title = "Scores plot of beef")
ggsave(paste(plot_path, "biplot_beef.png", sep = "/"), beef_biplot, device = "png", 
       height = 5, width = 7)


# loadings plots
fviz_pca_var(beef_pca, # factoextra
             geom = c("point"),
             col.var = "contrib",
             gradient.cols = viridis(n = 3, direction = -1),
             repel = TRUE)

new_beef_pca <- beef_pca
rownames(new_beef_pca$var$coord) <-short_lipid_names

plot_loadings(new_beef_pca, colour = TRUE, top_loadings = 10, xlab = "PC1", ylab = "PC2")
ggsave(paste(plot_path, "loadings1_beef.png", sep = "/"), 
       device = "png", 
       width = 7, 
       height =5)
plot_contrib_to_pc(new_beef_pca)
cpc1 <- fviz_contrib(new_beef_pca, choice = "var", axes = 1, top = 10,
                     fill = viridis(n = 1, begin = 0.3), color = viridis(n = 1, begin = 0.3),
                     title = "Contribution to PC1",
                     ggtheme = my_theme)
cpc2 <- fviz_contrib(new_beef_pca, choice = "var", axes = 2, top = 10,
                     fill = viridis(n = 1, begin = 0.3), color = viridis(n = 1, begin = 0.3),
                     title = "Contribution to PC2",
                     ggtheme = my_theme,
                     linecolor = "black", 
                     xtickslab.rt = 45)
contrib_pc <- ggarrange(plotlist = list(cpc1, cpc2), 
                        nrow = 1, 
                        ncol = 2, 
                        widths = c(1, 1), 
                        labels = c("A", "B"), 
                        font.label = list(size = 10, 
                                          color = "grey40", 
                                          face = "plain", 
                                          family = "AvantGarde"))
ggsave(paste(plot_path, "contrib_pc_beef.png", sep = "/"), 
       contrib_pc, 
       device = "png", 
       width = 10, 
       height = 5)



# most contributional lipids
meat_var <- beef_pca$var
meat_contrib <- as.data.frame(meat_var$contrib)
pc1_contrib_table <- meat_contrib[order(meat_contrib$Dim.1, decreasing = TRUE),]
pc1_contrib <- rownames(pc1_contrib_table)[1:10]
pc2_contrib_table <- meat_contrib[order(meat_contrib$Dim.2, decreasing = TRUE),]
pc2_contrib <- rownames(pc2_contrib_table)[1:10]
beef_pca_sub <- subset(lipid_list, select = c("Sample", "Group", pc1_contrib, pc2_contrib))

############### Clustering #################################################################
### Clustering
meat_clust <- data.frame(Group = lipid_list$group)
meat_clust <- cbind(meat_clust, select_if(lipid_list, is.numeric))
rownames(meat_clust) <- lipid_list$SID

ptab <- hclust_performance_table(meat_clust, 
                         dist_methods = c("euclidean", "maximum", "manhattan", "canberra", "binary",
                                          "minkowski"), 
                         hclust_methods = c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty",
                                            "median", "centroid"))
pplot <- hclust_performance_plot(meat_clust, 
                        dist_methods = c("euclidean", "maximum", "manhattan", "canberra",
                                         "minkowski"), 
                        hclust_methods = c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty",
                                           "median", "centroid"))

ggsave(paste(plot_path, "clust_performance.png", sep = "/"), 
       pplot, 
       device = "png", 
       width = 7, 
       height = 5)

meat_dist <- dist(select_if(lipid_list, is.numeric), method = "manhattan")
meat_hclust <- hclust(meat_dist, method = "average")
png(paste(plot_path, "dendrogram.png", sep = "/"), width = 20, height = 14, units = "cm", res = 3000)
hclust_dendrogram(meat_hclust,
                  labs = paste(lipid_list$sample,
                               meat_clust$Group, sep = "-"))
dev.off()

hclust_heatmap(lipid_list,
               dist_method = "manhattan",
               hclust_method = "average",
               row_names = meat_clust$Group, 
               out_path = paste(plot_path, "hclust_heat.png", sep = "/"))

new_lipid_list <- select_if(lipid_list, is.numeric)
new_lipid_list <- cbind(lipid_list$Group, new_lipid_list)

heatmap <- hclust_heatmap_interactive(new_lipid_list,
                           dist_method = "manhattan",
                           hclust_method = "average", 
                           row_names = lipid_list$Bio_replicate 
                           )

hclust_map <- hclust_heatmap_interactive(meat_clust,
                           dist_method = "manhattan",
                           hclust_method = "average", out_path = "heatmaply.png")


# t_lipid_list <- as.data.frame(t(select_if(lipid_list, is.numeric)))
# colnames(t_lipid_list) <- lipid_list$Bio_replicate
# 
# pplot2 <- hclust_performance_plot(t_lipid_list, 
#                                  dist_methods = c("euclidean", "maximum", "manhattan", "canberra",
#                                                   "minkowski"), 
#                                  hclust_methods = c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty",
#                                                     "median", "centroid"))
# 
# 
# lipid_hclust <- hclust(dist(t_lipid_list, method = "manhattan"), method = "average")
# hclust_dendrogram(lipid_hclust)

#################################################################################

############## Beef univariate ##################################################

p_beef <- one_sample_test_by_col(beef, beef$nutrition, method = t.test)
adj_beef <- p.adjust(p_beef$p_values, method = "fdr")
fc_beef <- log2_foldchange(beef,
                           beef$nutrition,
                           control_group = "stall",
                           test_group = "grazing")

beef_volcano <- data.frame(p_value = p_beef, adj_p_value = adj_beef, log2_foldchange = fc_beef)

vp <- volcano_plot(beef_volcano,
             foldchange_col = beef_volcano$log2_foldchange,
             significance_col = beef_volcano$adj_p_value,
             foldchange = 1,
             significance = 0.05)

ggsave(filename = paste(plot_path, "beef_volcano.png", sep = "/"), 
       plot = vp, 
       device = "png", 
       width = 10, 
       height = 5
)
