### Install packages
# install.packages("ggplot2")
# install.packages("stringr")
# install.packages("tidyr")
# install.packages("data.table")
# install.packages("textshape")
## BiocManager::install("mixOmics")
## BiocManager::install("RVAideMemoire") 
# install.packages("MASS")
# install.packages("psych")
# install.packages("dplyr")
# install.packages("devtools")
# devtools::install_github("ricardo-bion/ggradar", dependencies=TRUE)
# install.packages("fmsb")
# install.packages("RColorBrewer")
# install.packages("scales")


### load packages
# library(gridExtra)
library(stringr) # count separators
library(ggplot2)#, # plots
#library(ggradar, scales) # radar chart with ggplot
# library(tidyr)
library(data.table) # transpose data frame
# library(textshape)
library(tibble) # data frame manipulation
# library(RVAideMemoire) 
library(MASS) # for paralell plot 
library(psych) # for correlation plot 
library(dplyr) # select part of data
library(fmsb) # spider chart
library(RColorBrewer) # pretty color combinations
library(scales) # scale opacity of filling (alpha)

source("lipidome_comparison_functions.R")

# set ggplot theme
theme_set(
  theme_minimal() +
    theme(legend.position = "top")
)

## set variables
working_directory <- "/home/lisa/FH/Masterarbeit/LipidomeComparison"
setwd(working_directory)

input_path <- "/home/lisa/FH/Masterarbeit/LipidomeComparison/data/Probe-Datensatz_lisa.csv"
test_path <- "/home/lisa/FH/Masterarbeit/LipidomeComparison/data/test2.csv"
plot_path <- paste(working_directory, "/plots", sep = "")
plot_name <- paste(plot_path, "/test_data", sep = "")

## load & transform data
lipid_data <- read.csv(input_path, sep = ",", dec = ".", header = TRUE) #read data
test_data <- read.csv(test_path, sep = ",", dec = ".", header = TRUE)

t_lipid_data <- pretty_transpose(lipid_data)
t_test_data <- pretty_transpose(test_data)

working_data <- t_lipid_data
working_data <- SID_to_metadata(working_data)
working_data <- character_to_factor(working_data) 


## summary biological & technical replicates
means_biol <- calc_by_replicate(working_data, "treatment", mean)
means_tech <- calc_by_replicate(working_data, "biol_replicate", mean)

sd_biol <- calc_by_replicate(working_data, "treatment", sd)
sd_tech <- calc_by_replicate(working_data, "biol_replicate", sd)

## plots for normal distribution
qqplot_by_factor(working_data, "treatment", plot_name)
histogram_by_factor(working_data, "treatment", plot_name)
boxplot_by_factor(working_data, "treatment", plot_name)

## test for normal distribution
### Don't use with multi modal data --> check histogram and qq plots first
shapiro_all <- lapply(dplyr::select_if(working_data, is.numeric), shapiro.test)
shapiro_all <- sapply(shapiro_all, `[`, c("statistic","p.value"))

shapiro_by_treatment <- shapiro_by_factor(working_data, "treatment")

## check for correlations between lipids
correlation_plot(working_data, "pearson") # for <= 10 variables

### plots 
## paralell plot for <= 10 variables
parallel_plot(working_data, "biol_replicate", plot_name) 

## spider chart
spider_data <- SID_to_metadata(t_lipid_data) # calculate means so there is only one value per group
spider_data <- calc_by_replicate(spider_data, "treatment", mean)

spider_chart(spider_data)


### PCA
lipid_pca <- prcomp(select_if(working_data, is.numeric), scale = FALSE, center = TRUE)
summary(lipid_pca)

{plot(lipid_pca,
      main = NULL)
  title(main = NULL)}

biplot(lipid_pca)
