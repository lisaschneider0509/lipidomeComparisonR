### load packages
library(data.table) # transpose data frame
library(tidyverse)
library(viridis) # colorblind save color schemes
library(ggpubr) # multiple plots on one page
library(ggmosaic)

source("R/lipidome_comparison_dataTransformaions.R")

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
  
  working_directory <- "/home/lisa/FH/Masterarbeit/LipidomeComparison"
  setwd(working_directory)
  
  data_dir <- "/home/lisa/FH/Masterarbeit/LipidomeComparison/data"
  lipid_list_path <- "/home/lisa/FH/Masterarbeit/LipidomeComparison/data/meat_fish_final_raw.csv"
  annotation_path <- "/home/lisa/FH/Masterarbeit/LipidomeComparison/data/meat_annotation.csv"
  data_matrix_path <- paste("/home/lisa/FH/Masterarbeit/LipidomeComparison/data/", project, "_data_matrix.csv", sep = "")
  
  plot_path <- paste(working_directory, "/plots", sep = "")

########### Read imputed data #############
  complete_data <- read.csv(paste(data_dir, "/", project, "_complete_data.csv", sep = ""), row.names = 1)
  colnames(complete_data) <- unlist(read.csv(paste(data_dir, "/", project, "_complete_data.csv", sep = ""), header = FALSE)[1, -1], use.names = FALSE)
  
  lipids <- colnames(dplyr::select_if(complete_data, is.numeric))
  
########### Get lipid class ###############
  lipid_class <- sub("\\ .*","", lipids)
  lipid_class <- sub(".*\\-", "", lipid_class)
  
########### Get lipid species #############
  no_class <- sub(".*\\(", "", lipids)
  lipid_species <- sub("\\-.*", "", no_class)
  lipid_species <- sub("\\).*", "", lipid_species)
  
########### Get fatty acyl chains #########
  fatty_acyls <- sub(".*\\-", "", no_class)
  fatty_acyls <- sub("\\).*", "", fatty_acyls)
  fatty_acyls <- sub(".*e.*", NA, fatty_acyls)
  fatty_acyls <- ifelse(is.na(fatty_acyls), lipid_species, fatty_acyls)
  
  fatty_acyls <- sub("P", "P-", fatty_acyls)
  fatty_acyls <- sub("O", "O-", fatty_acyls)
  
########## Get ether bonds ################
  alkenyl_ether <- grepl("P|O", fatty_acyls)
  
########## Get total chain length #########
  total_chain_length <- as.numeric(sub("\\:.*", "", lipid_species))
  
########## Get total number of unsaturated bonds ###
  unsaturated_bonds <- as.numeric(sub(".*\\:", "", lipid_species))
  
  is_pufa <- unsaturated_bonds >= 2 # check if lipid contains >= 1 PUFA (chek for two MUFA manually)
  
########## Data frame with lipid information ###
  
  lipid_information <- data.frame(lipidID = lipids, 
                           lipidClass = lipid_class, 
                           totalChainLength = total_chain_length,
                           unsaturation = unsaturated_bonds,
                           lipidSpecies = paste(lipid_class, lipid_species), 
                           fattyAcyls = paste(lipid_class, fatty_acyls), 
                           PUFA = is_pufa, 
                           ether = alkenyl_ether)
  
######### Plot lipid information #########
  
  class_table <- as.data.frame(table(lipid_class))
  species_table <- as.data.frame(table(sub(".*\\ ", "", lipid_information$lipidSpecies)))
  alcyl_table <- as.data.frame(table(sub(".*\\ ", "", lipid_information$fattyAcyls)))
  
  class_plot <- ggplot(data = lipid_information) +
    geom_mosaic(aes(x = product(lipidClass), fill=lipidClass)) +
    scale_fill_viridis_d(alpha = 1) +
    labs(x = "lipid class", y = NULL, title='Lipid classes') +
    theme(axis.text.x=element_text(angle=90), #
          axis.text.y = element_blank(),
          legend.position = "bottom", 
          legend.title = element_blank())

  pufa_plot <- ggplot(data = lipid_information) +
    geom_mosaic(aes(x = product(PUFA), fill=PUFA)) +
    scale_fill_viridis_d(name = "", labels = c("no PUFA", "PUFA")) +
    labs(x=NULL, y = NULL, title='Lipids with at least one PUFA') + 
    theme(legend.position = "bottom", 
          axis.text.y = element_blank(),
          legend.title = element_blank())
  
  freq_plot <- ggpubr::ggarrange(plotlist = list(class_plot, pufa_plot),
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
  
####### Rename data in LIPID MAPS nomenclature ####
  short_lipid_names <- paste(lipid_information$fattyAcyls, " (", rowid(lipid_information$fattyAcyls), ")", sep = "")
  colnames(complete_data) <- c(colnames(complete_data)[1], 
                               short_lipid_names, 
                               colnames(complete_data)[(ncol(complete_data)-1):ncol(complete_data)])

write.csv(file = paste(data_dir, "/", project, "_renamed_data.csv", sep = ""), x = complete_data)
write.csv(file = paste(data_dir, "/", project, "_lipid_information.csv", sep = ""), x = lipid_information)

