### load packages
{library(dplyr) # select part of data
  library(stringr) # count separators
  library(data.table) # transpose data frame
  # library(crmn)
  # library(impute)
  # library(imputeLCMD)
  library(tibble) # data frame manipulation
  
  library(tidyverse)
  library(ggplot2)#, # plots
  library(viridis) # colorblind save color schemes
  # library(GGally) # paralell plot
  # library(fmsb) # spider chart
  # library(scales) # scale opacity of filling (alpha)
  # library(ggpubr) # multiple plots on one page
  # library(ggmosaic)
  
  # library(ggrepel)
  # library(factoextra)
  # library(ggfortify) # biplot with ggplot
  # library(corrplot)
  # library(FactoMineR)
  
  library(heatmaply) # interactive heatmap
  library(gplots) # heatmap
  library(plotly) # interactive ggplots
  library(htmlwidgets) # save plotly-plots as html
  library(dendextend)
  
  source("R/lipidome_comparison_dataTransformaions.R")
  source("R/lipidome_comparison_EDA.R")
  source("R/lipidome_comparison_clustering.R")
  source("R/lipidome_comparison_hypothesis_testing.R")}

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

############## import lipid data #############################

lipid_data <- read.csv(file = paste(data_dir, "/", project, "_renamed_data.csv", sep = ""), row.names = 1)
colnames(lipid_data) <- base::unlist(read.csv(paste(data_dir, "/", project, "_renamed_data.csv", sep = ""), header = FALSE)[1, -1], use.names = FALSE)
############### Clustering #################################################################
### Clustering
meat_clust <- data.frame(Group = lipid_data$Group)
meat_clust <- cbind(meat_clust, select_if(lipid_data, is.numeric))
rownames(meat_clust) <- lipid_data$Replicate

ptab <- hclust_performance_table(meat_clust, 
                                 dist_methods = c("euclidean", "maximum", "manhattan", "canberra", "minkowski"), 
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
       width = 10, 
       height = 5)

meat_dist <- dist(select_if(lipid_data, is.numeric), method = "manhattan")
meat_hclust <- hclust(meat_dist, method = "average")
png(paste(plot_path, "dendrogram.png", sep = "/"), width = 20, height = 20, units = "cm", res = 200)
hclust_dendrogram(meat_hclust,
                  labs = paste(lipid_data$Sample,
                               meat_clust$Group, sep = "-"))
dev.off()

hclust_heatmap(lipid_data,
               dist_method = "manhattan",
               hclust_method = "average",
               row_names = meat_clust$Group, 
               out_path = paste(plot_path, "hclust_heat.png", sep = "/"))

new_lipid_data <- select_if(lipid_data, is.numeric)
new_lipid_data <- cbind(lipid_data$Group, new_lipid_data)

heatmap <- hclust_heatmap_interactive(new_lipid_data,
                                      dist_method = "manhattan",
                                      hclust_method = "average", 
                                      row_names = lipid_data$Bio_replicate 
)

hclust_map <- hclust_heatmap_interactive(meat_clust,
                                         dist_method = "manhattan",
                                         hclust_method = "average", out_path = "heatmaply.png")


# t_lipid_data <- as.data.frame(t(select_if(lipid_data, is.numeric)))
# colnames(t_lipid_data) <- lipid_data$Bio_replicate
# 
# pplot2 <- hclust_performance_plot(t_lipid_data, 
#                                  dist_methods = c("euclidean", "maximum", "manhattan", "canberra",
#                                                   "minkowski"), 
#                                  hclust_methods = c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty",
#                                                     "median", "centroid"))
# 
# 
# lipid_hclust <- hclust(dist(t_lipid_data, method = "manhattan"), method = "average")
# hclust_dendrogram(lipid_hclust)

