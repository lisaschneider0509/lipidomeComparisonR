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
  
  source("R/lipidome_comparison_dataTransformaions.R")
  source("R/lipidome_comparison_EDA.R")
  source("R/lipidome_comparison_pca.R")
  source("R/lipidome_comparison_clustering.R")
  source("R/lipidome_comparison_hypothesis_testing.R")
  }

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

annotation_data <- read.csv(annotation_path, sep = ",", dec = ".", header = TRUE)

map <- data.frame(Sample=annotation_data$Sample,
                  Nurtrition=annotation_data$Nutrition)
lipid_data <- left_join(lipid_data, map, by="Sample")

beef <- subset(lipid_data, subset = lipid_data$Group == "beef")
beef <- droplevels(beef)

############## Beef univariate ##################################################

p_beef <- one_sample_test_by_col(beef, beef$Nurtrition, method = t.test)
adj_beef <- p.adjust(p_beef$p_values, method = "fdr")
fc_beef <- log2_foldchange(beef,
                           beef$Nutrition,
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
