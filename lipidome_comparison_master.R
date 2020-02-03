### Install packages
# install.packages("ggplot2")
# install.packages("stringr")
# install.packages("tidyr")
# install.packages("data.table")
# install.packages("textshape")
BiocManager::install("mixOmics")
BiocManager::install("RVAideMemoire")

### load packages
library(gridExtra)
library(stringr)
library(ggplot2)
library(tidyr)
library(data.table)
library(textshape)
library(tibble)
library(RVAideMemoire)

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
t_lipid_data <- read_transpose(input_path)
t_test_data <- read_transpose(test_path)

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
# shapiro_all <- lapply(working_data[-(1:3)], shapiro.test)
# shapiro_all <- sapply(shapiro_all, `[`, c("statistic","p.value"))

aggregate(working_data, by = list(working_data[[2]]),
          FUN = function(x) {y <- shapiro.test(x); c(y$statistic, y$p.value)})

# RVAideMemoire::byf.shapiro(working_data$`5-HETE` ~ working_data$treatment, working_data)
