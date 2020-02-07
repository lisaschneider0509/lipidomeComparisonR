### Install packages
## general
# install.packages("tibble")
# install.packages("stringr")
# install.packages("data.table")
# install.packages("dplyr")
# install.packages("devtools")

## graphs
# install.packages("ggplot2")
# install.packages("scales")
# install.packages("viridis")

## correlation plot 
# install.packages("psych")

## for spider chart
# install.packages("fmsb")

## for paralell plot
# install.packages("GGally")

# install.packages("hrbrthemes")

## for PCA
# install.packages("ggfortify")



### load packages
# library(gridExtra)
library(stringr) # count separators
library(ggplot2)#, # plots
library(data.table) # transpose data frame
library(tibble) # data frame manipulation
library(viridis) # colorblind save color schemes
library(GGally, hrbrthemes) # paralell plot
library(psych) # for correlation plot 
library(dplyr) # select part of data
library(fmsb) # spider chart
library(scales) # scale opacity of filling (alpha)
library(devtools)
library(ggfortify)

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
parallel_plot(working_data[1:10], "treatment", plot_name) 

## spider chart
spider_data <- SID_to_metadata(t_lipid_data) # calculate means so there is only one value per group
spider_data <- calc_by_replicate(spider_data, "treatment", mean)
rownames(spider_data) <- spider_data$Group.1

spider_chart(spider_data[1:10])


### PCA
wd <- working_data
lipid_pca <- prcomp(select_if(wd[1:30], is.numeric), scale = FALSE, center = TRUE)
summary(lipid_pca)
# 
# {plot(lipid_pca,
#       main = NULL)
#   title(main = NULL)}
# 
# biplot(lipid_pca)

autoplot(lipid_pca, data = wd, colour = 'treatment',
         loadings = TRUE,
         loadings.label = TRUE, 
         loadings.label.size = 3, 
         frame = TRUE, 
         frame.type = "norm"
         ) + scale_color_viridis(discrete=TRUE)


