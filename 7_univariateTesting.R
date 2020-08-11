if (!requireNamespace("tidyverse", quietly = TRUE)){
  install.packages("tidyverse")}
if (!requireNamespace("ggpubr", quietly = TRUE)){
  install.packages("ggpubr")}

if (!requireNamespace("lipidomeComparison", quietly = TRUE)){
  devtools::install_local("lipidomeComparison_0.1.0.tar.gz")}  


### load packages
library(tidyverse)
library(ggpubr)
library(lipidomeComparison)


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
project <- "meat"

working_directory <- "/home/lisa/FH/Masterarbeit/meatLipidomics"
setwd(working_directory)

data_dir <- "/home/lisa/FH/Masterarbeit/meatLipidomics/data"
lipid_list_path <- "/home/lisa/FH/Masterarbeit/meatLipidomics/data/meat_fish_final_raw.csv"
annotation_path <- "/home/lisa/FH/Masterarbeit/meatLipidomics/data/meat_annotation.csv"
data_matrix_path <- paste("/home/lisa/FH/Masterarbeit/meatLipidomics/data/", project, "_data_matrix.csv", sep = "")

plot_path <- paste(working_directory, "/plots", sep = "")


############## import lipid data & subset to two groups #############################

lipid_data <- read.csv(file = paste(data_dir, "/", project, "_renamed_data.csv", sep = ""), row.names = 1)
colnames(lipid_data) <- base::unlist(read.csv(paste(data_dir, "/", project, "_renamed_data.csv", sep = ""), header = FALSE)[1, -1], use.names = FALSE)

annotation_data <- read.csv(annotation_path, sep = ",", dec = ".", header = TRUE)

map <- data.frame(Sample=annotation_data$Sample,
                  Nutrition=annotation_data$Nutrition)
lipid_data <- left_join(lipid_data, map, by="Sample")

beef <- subset(lipid_data, subset = lipid_data$Group == "beef")
beef <- droplevels(beef)

############## Beef univariate ##################################################

p_beef <- one_sample_test_by_col(beef, beef$Nutrition, method = t.test)
adj_beef <- p.adjust(p_beef$p_values, method = "fdr")
fc_beef <- log2_foldchange(beef,
                           beef$Nutrition,
                           control_group = "stall",
                           test_group = "grazing")

beef_volcano <- data.frame(p_value = p_beef, adj_p_value = adj_beef, log2_foldchange = fc_beef)

volcano_plot <- volcano_plot(beef_volcano,
                   foldchange_col = beef_volcano$log2_foldchange,
                   significance_col = beef_volcano$adj_p_value,
                   foldchange = 1,
                   significance = 0.05)

ggsave(filename = paste(plot_path, "beef_volcano.png", sep = "/"), 
       plot = volcano_plot, 
       device = "png", 
       width = 10, 
       height = 5
)
