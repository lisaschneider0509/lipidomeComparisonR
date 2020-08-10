### load packages
{library(dplyr) # select part of data
  library(stringr) # count separators
  library(data.table) # transpose data frame
  library(tibble) # data frame manipulation
  
  library(tidyverse)
  library(ggplot2)#, # plots
  library(viridis) # colorblind save color schemes
  library(GGally) # paralell plot
  library(fmsb) # spider chart
  library(scales) # scale opacity of filling (alpha)
  library(ggpubr) # multiple plots on one page
  # library(ggmosaic)
  
  library(ggrepel)
  library(factoextra)
  library(ggfortify) # biplot with ggplot
  library(FactoMineR)
  
  source("R/lipidome_comparison_dataTransformaions.R")
  source("R/lipidome_comparison_EDA.R")
  source("R/lipidome_comparison_pca.R")
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

###################### PCA #############################################
  lipid_pca <- FactoMineR::PCA(select_if(lipid_data, is.numeric), scale.unit = T, graph = F)
  
############ get eigenvalues and plot contribution of PCs to variance #####
  lipid_eigenvalues <- as.data.frame(factoextra::get_eigenvalue(lipid_pca))

{ # scree
  # eigen_plot <- plot_pc_variance(lipid_eigenvalues[1:10,], 
  #                                x = seq(1:10), 
  #                                y = lipid_eigenvalues$eigenvalue[1:10], 
  #                                title = "Eigenvalue", 
  #                                ylab = "eigenvalues")
  # 
  var_plot <- plot_pc_variance(lipid_eigenvalues[1:10,], 
                               x = seq(1:10), 
                               y = lipid_eigenvalues$variance.percent[1:10], 
                               xlab = "principal components",
                               title = "Variance [%]")
  
  cum_var_plot <- plot_pc_variance(lipid_eigenvalues[1:10,], 
                                   x = seq(1:10), 
                                   y = lipid_eigenvalues$cumulative.variance.percent[1:10], 
                                   title = "Cummulative variance [%]", 
                                   ylab = NULL,
                                   xlab = "principal components",
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
  ggsave(paste(plot_path, "/", project, "_scree.png", sep = ""), scree, device = "png", 
         height = 5, width = 10)
}

################# plot sample scores ###############################
lipid_biplot <- biplot_ggplot2(lipid_data, 
                               groups = "Group", 
                               loadings = FALSE, 
                               ellipse = TRUE, 
                               scale = TRUE, 
                               title = "Scores plot of meat data")
ggsave(filename = paste(plot_path, "/", project, "_biplot.png", sep = ""), 
       plot = lipid_biplot, 
       device = "png", 
       height = 5, width = 10)


################# plot variable loadings ##########################
fviz_pca_var(lipid_pca, # factoextra
             geom = c("point"),
             col.var = "contrib",
             gradient.cols = viridis(n = 3, direction = -1),
             repel = TRUE)


loadings_plot <- plot_loadings(lipid_pca, 
                               colour = TRUE, 
                               top_loadings = 10, 
                               xlab = "PC1", 
                               ylab = "PC2", 
                               title = "Loadings of meat data")

ggsave(paste(plot_path, "/", project, "_loadings_plot.png", sep = ""), 
       plot = loadings_plot,
       device = "png", 
       width = 10, 
       height = 5)

################# plot contribution of variables to principal components ############

contribution_to_pc <- plot_contrib_to_pc(lipid_pca)

pc1_loadngs_bar <- fviz_contrib(lipid_pca, choice = "var", axes = 1, top = 10,
                     fill = viridis(n = 1, begin = 0.3), color = viridis(n = 1, begin = 0.3),
                     title = "Contribution to PC1",
                     ggtheme = my_theme)

pc2_loadings_bar <- fviz_contrib(lipid_pca, choice = "var", axes = 2, top = 10,
                     fill = viridis(n = 1, begin = 0.3), color = viridis(n = 1, begin = 0.3),
                     title = "Contribution to PC2",
                     ggtheme = my_theme,
                     linecolor = "black", 
                     xtickslab.rt = 45)

contribution_plot <- ggarrange(plotlist = list(pc1_loadngs_bar, pc2_loadings_bar), 
                        nrow = 1, 
                        ncol = 2, 
                        widths = c(1, 1), 
                        labels = c("A", "B"), 
                        font.label = list(size = 10, 
                                          color = "grey40", 
                                          face = "plain", 
                                          family = "AvantGarde"))

ggsave(paste(plot_path, "/", project, "_pc_contribution.png", sep = ""), 
       contribution_plot, 
       device = "png", 
       width = 10, 
       height = 5)


###################### get the most contributional variables ##################
lipid_var <- lipid_pca$var
lipid_contrib <- as.data.frame(lipid_var$contrib)
pc1_contrib_table <- lipid_contrib[order(lipid_contrib$Dim.1, decreasing = TRUE),]
pc1_contrib <- rownames(pc1_contrib_table)[1:10]
pc2_contrib_table <- lipid_contrib[order(lipid_contrib$Dim.2, decreasing = TRUE),]
pc2_contrib <- rownames(pc2_contrib_table)[1:10]
lipid_pca_sub <- subset(lipid_data, select = c("Replicate", "Group", pc1_contrib, pc2_contrib))

pp <- parallel_plot(lipid_pca_sub, lipid_pca_sub$Group, title = "Relative abundances of most contributional lipids")
ggsave(filename = paste(plot_path, "/", project, "_paralellPlot.png", sep = ""), plot = pp, 
       width = 10, 
       height = 5)

