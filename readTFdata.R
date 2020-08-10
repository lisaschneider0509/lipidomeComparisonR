### load packages
{library(dplyr) # select part of data
  library(stringr) # count separators
  library(data.table) # transpose data frame
  library(crmn)
  library(impute)
  library(imputeLCMD)
  library(tibble) # data frame manipulation
  
  # library(tidyverse)
  # library(ggplot2)#, # plots
  # library(viridis) # colorblind save color schemes
  # library(GGally) # paralell plot
  # library(fmsb) # spider chart
  # library(scales) # scale opacity of filling (alpha)
  # library(ggpubr) # multiple plots on one page
  # library(ggmosaic)
  # 
  source("R/lipidome_comparison_dataTransformaions.R")
  source("R/lipidome_comparison_EDA.R")}

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
                    Group=annotation_data$Meat.type)
  target_lipids <- left_join(target_lipids, map, by="Sample_nr")
  # target_lipids <- target_lipids[-ncol(target_lipids)]
  

write.csv(x = target_lipids, paste(data_dir, "/", project, "_condensed_data.csv", sep = ""))

