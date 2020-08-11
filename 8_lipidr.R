### install packages
if (!requireNamespace("tidyverse", quietly = TRUE)){
  install.packages("tidyverse")}
if (!requireNamespace("ggpubr", quietly = TRUE)){
  install.packages("ggpubr")}
if (!requireNamespace("viridis", quietly = TRUE)){
  install.packages("viridis")}
if (!requireNamespace("BiocManager::", quietly = TRUE)){
  install.packages("BiocManager::")}
if (!requireNamespace("lipidr", quietly = TRUE)){
  BiocManager::install("lipidr")}

if (!requireNamespace("lipidomeComparison", quietly = TRUE)){
  devtools::install_local("lipidomeComparison_0.1.0.tar.gz")} 


library(tidyverse)
library(lipidr)
library(ggpubr)
library(viridis)


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



################### read data ###############
lipidr_list <- read.csv(file = paste(data_dir, "/", project, "_renamed_data.csv", sep = ""), row.names = 1)
colnames(lipidr_list) <- base::unlist(read.csv(paste(data_dir, "/", project, "_renamed_data.csv", sep = ""), header = FALSE)[1, -1], use.names = FALSE)
lipidr_list <- subset(lipidr_list, subset = lipidr_list$Group == "beef")


########### transform to lipidomics experiment data type ################
lipidr_beef <- subset(lipidr_list, lipidr_list$Group == "beef")
droplevels(lipidr_beef)

lipidr_matrix <- t(select_if(lipidr_beef, is.numeric))
colnames(lipidr_matrix) <- lipidr_beef$Replicate

annotation_df <- subset(lipidr_beef, select = c("Replicate", "Sample"))
annotation_data <- read.csv(annotation_path, sep = ",", dec = ".", header = TRUE)
map <- data.frame(Sample = annotation_data$Sample, 
                  Nutrition = annotation_data$Nutrition)
annotation_df <- left_join(annotation_df, map, by = "Sample")
names(annotation_df) <- c("Sample", "ID", "Nutrition")

lipidr_experiment <- as_lipidomics_experiment(lipidr_matrix)
lipidr_experiment <- add_sample_annotation(lipidr_experiment, annotation_df)


############ perform and visualize univariate analysis ##############
de_results = de_analysis(
  data=lipidr_experiment, 
  group_col = "Nutrition",
  stall - grazing,
  measure="Area")
head(de_results)


vp <- volcano_plot(volcano_df = de_results, foldchange_col = de_results$logFC, significance_col = de_results$adj.P.Val)
ggsave(filename = paste(plot_path, "beef_volcano.png", sep = "/"), 
       plot = vp, 
       device = "png", 
       width = 7, 
       height = 5
)

subset(de_results$Molecule, de_results$adj.P.Val <= 0.05)

############## lipid set enrichment analysis #####################
enrich_results <- lipidr::lsea(de_results, rank.by = "logFC")
sig_lipids <- lipidr::significant_lipidsets(enrich_results)
sig_lipids <- sub(".*_", "", sig_lipids$`stall - grazing`)

lipidClass <- ggplot(data = de_results) + 
  geom_boxplot(aes(x = Class, y = logFC, fill = is.element(de_results$Class, sig_lipids))) +
  scale_fill_manual(name = 'Significant', values = setNames(c(viridis(n = 1, begin = 0.5),'grey90'),c(T, F))) +
  labs(title = "Foldchange by lipid classes", y = "log2FC", x = "lipid class")

totalChain <- ggplot(data = de_results) + 
  geom_boxplot(aes(x = as.factor(total_cl), y = logFC, fill = is.element(de_results$total_cl, sig_lipids))) +
  scale_fill_manual(name = 'Significant', values = setNames(c(viridis(n = 1, begin = 0.5),'grey90'),c(T, F))) +
  labs(title = "Foldchange by total chain length", y = "log2FC", x = "total chain length")

unsaturatedBonds <- ggplot(data = de_results) + 
  geom_boxplot(aes(x = as.factor(total_cs), y = logFC, fill = is.element(de_results$total_cs, sig_lipids))) +
  scale_fill_manual(name = 'Significant', values = setNames(c(viridis(n = 1, begin = 0.5),'grey90'),c(T, F))) +
  labs(title = "Foldchange by unsaturated bonds", y = "log2FC", x = "total unsaturated bonds")


mylist <- list(lipidClass, totalChain, unsaturatedBonds)

lsea <- ggarrange(plotlist = mylist, ncol = 3, nrow = 1, common.legend = TRUE, legend = "bottom", 
                  labels = c("(A)", "(B)", "(C)"), 
                  font.label = list(size = 10, 
                                    color = "grey40", 
                                    face = "plain", 
                                    family = "AvantGarde"))
ggsave(filename = paste(plot_path, "lsea.png", sep = "/"), plot = lsea, device = "png", 
       width = 15, height = 5)
