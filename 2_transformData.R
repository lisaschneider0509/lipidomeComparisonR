### load packages
library(data.table) # transpose data frame
library(impute)
library(imputeLCMD)
library(tibble) # data frame manipulation
library(tidyverse)

source("R/lipidome_comparison_dataTransformaions.R")
  # source("R/lipidome_comparison_EDA.R")

############### set variables #############################
project <- "meat"

working_directory <- "/home/lisa/FH/Masterarbeit/LipidomeComparison"
setwd(working_directory)

data_dir <- "/home/lisa/FH/Masterarbeit/LipidomeComparison/data"
lipid_list_path <- "/home/lisa/FH/Masterarbeit/LipidomeComparison/data/meat_fish_final_raw.csv"
annotation_path <- "/home/lisa/FH/Masterarbeit/LipidomeComparison/data/meat_annotation.csv"

plot_path <- paste(working_directory, "/plots", sep = "")

######### Load condensed data #####################

condensed_data <- read.csv(paste(data_dir, "/", project, "_condensed_data.csv", sep = ""), row.names = 1)
colnames(condensed_data) <- unlist(read.csv(paste(data_dir, "/", project, "_condensed_data.csv", sep = ""), header = FALSE)[1, -1], use.names = FALSE)
  
annotation_data <- read.csv(annotation_path)
  
######### Log 2 transformation ####################
log2_data <- log2(dplyr::select_if(condensed_data, is.numeric))
log2_data <- cbind(dplyr::select_if(condensed_data, is.character), 
                   dplyr::select_if(condensed_data, is.factor), 
                   log2_data)

######### Remove columns with >20% missing values ###
  log2_data <- log2_data[, which(colMeans(!is.na(log2_data)) > 0.8)]
  log2_data <- as.matrix(dplyr::select_if(log2_data, is.numeric))
  
######### Data imputation ###########################
  imputed_data <- as.data.frame(impute.QRILC(log2_data, tune.sigma = 1)[[1]])

  
######### Normaization ##############################
  # normalized_data <- med_normalize(imputed_data)
  # numeric_data <- normalized_data
  
  ### no normalization 
  numeric_data <- imputed_data
  
  #### calculate the means for the replicates
  data_groups <- generate_categorical_table(condensed_data$Group)
  complete_data <- calc_by_replicate(numeric_data, condensed_data$Biol_rep, mean)
  colnames(complete_data)[1] <- "Replicate"
  complete_data$Sample <- sub("\\_.*", "", complete_data$Replicate)
  
  map <- data.frame(Sample=annotation_data$Sample,
                    Group=annotation_data$Meat.type)
  complete_data <- dplyr::left_join(complete_data, map, by="Sample")

write.csv(file = paste(data_dir, "/", project, "_complete_data.csv", sep = ""), x = complete_data)
