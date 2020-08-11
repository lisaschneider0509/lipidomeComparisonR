### install packages
if (!requireNamespace("tidyverse", quietly = TRUE)){
  install.packages("tidyverse")}
if (!requireNamespace("data.table", quietly = TRUE)){
  install.packages("data.table")}
if (!requireNamespace("lipidomeComparison", quietly = TRUE)){
  devtools::install_local("lipidomeComparison_0.1.0.tar.gz")}  

### load packages
library(tidyverse)
library(data.table)
library(lipidomeComparison)
  
  

############### set variables #############################
project <- "meat"

working_directory <- "/home/lisa/FH/Masterarbeit/meatLipidomics"

data_dir <- "/home/lisa/FH/Masterarbeit/meatLipidomics/data"
lipid_list_path <- "/home/lisa/FH/Masterarbeit/meatLipidomics/data/meat_fish_final_raw.csv"
annotation_path <- "/home/lisa/FH/Masterarbeit/meatLipidomics/data/meat_annotation.csv"
data_matrix_path <- paste("/home/lisa/FH/Masterarbeit/meatLipidomics/data/", project, "_data_matrix.csv", sep = "")

plot_path <- paste(working_directory, "/plots", sep = "")


############## import TF data #############################

  lipid_list <- read.csv(lipid_list_path, sep = ",", dec = ".", header = TRUE)
  annotation_data <- read.csv(annotation_path, sep = ",", dec = ".", header = TRUE)
  
  lipid_list <- subset(lipid_list, select = c(Compound, RT, Type, Filename, Status, Area))
  lipid_list <- subset(lipid_list, Status == "Processed")
  lipid_list[lipid_list=='N/F'] <- NA
  lipid_list$Area <- as.numeric(as.character(lipid_list$Area))

  target_lipids <- subset(lipid_list, Type == "Target Compound")
  lipid_standards <- subset(lipid_list, Type == "Internal Standard")
  
  target_lipids <- flip_df(target_lipids)
  lipid_standards <- flip_df(lipid_standards)

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


write.csv(x = target_lipids, paste(data_dir, "/", project, "_condensed_data.csv", sep = ""))

