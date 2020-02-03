### Install packages
# install.packages("ggplot2")
# install.packages("stringr")
# install.packages("tidyr")
# install.packages("data.table")
# install.packages("textshape")
# BiocManager::install("mixOmics")
# BiocManager::install("RVAideMemoire")
# install.packages("MASS")
# install.packages("psych")
# install.packages("dplyr")

### load packages
library(gridExtra)
library(stringr)
library(ggplot2)
library(tidyr)
library(data.table)
library(textshape)
library(tibble)
library(RVAideMemoire)
library(MASS)
library(psych)
library(dplyr)

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

t_lipid_data <- read_transpose(lipid_data)
t_test_data <- read_transpose(test_data)

working_data <- t_test_data

working_data <- SID_to_metadata(working_data)



## summary biological & technical replicates

means_biol <- calc_biol_rep(working_data, mean)
means_tech <- calc_tech_rep(working_data, mean)

sd_biol <- calc_biol_rep(working_data, sd)
sd_tech <- calc_tech_rep(working_data, sd)


## plots for normal distribution
qqplot_by_factor(working_data, "treatment", plot_name)
histogram_by_factor(working_data, "treatment", plot_name)
boxplot_by_factor(working_data, "treatment", plot_name)


## test for normal distribution
### Don't use with multi modal data --> check histogram and qq plots first
shapiro_all <- lapply(working_data_numeric, shapiro.test)
shapiro_all <- sapply(shapiro_all, `[`, c("statistic","p.value"))

shapiro_by_treatment <- shapiro_by_factor(working_data, "treatment")


## check for correlations between lipids
correlation_plot(working_data, "pearson") # for <= 10 variables


### plots 
## paralell plot 
mycolors <- colors()[as.numeric(working_data$treatment)*11]
MASS::parcoord(dplyr::select_if(working_data, is.numeric), col = mycolors)

