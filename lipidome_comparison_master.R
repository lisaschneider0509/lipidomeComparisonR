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
library(factoextra)
library(plotly) # interactive heatmap
library(heatmaply) # interactive heatmap
library(gplots) # heatmap 
library(dendextend)
library(limma) # hypothesis testing

source("lipidome_comparison_dataTransformaions.R")
source("lipidome_comparison_EDA.R")
source("lipidome_comparison_pca.R")
source("lipidome_comparison_clustering.R")


# set ggplot theme
my_theme <- theme_set(
  theme_minimal() +
    theme(plot.title = element_text(size=12, hjust = 0.5, family="AvantGarde"),
          axis.text.x = element_text(size = 8, colour = "grey40", family="AvantGarde"),
          axis.text.y = element_text(size = 8, colour = "grey40", family="AvantGarde"),
          # axis.title = element_text(size = 10, colour = "grey40", family="AvantGarde"),
          axis.title.x = element_text(size = 10, hjust = 0.5, colour = "grey40", family="AvantGarde"),
          axis.title.y = element_text(size = 10, hjust = 0.5, colour = "grey40", family="AvantGarde"),
          legend.text = element_text(size = 8, colour = "grey40", family="AvantGarde"),
          legend.title = element_text(size = 10, colour = "grey40", family="AvantGarde"))
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
colnames(t_lipid_data) <- gsub(" ", "_", colnames(t_lipid_data), fixed = TRUE)
t_test_data <- pretty_transpose(test_data)
colnames(t_test_data) <- gsub(" ", "_", colnames(t_test_data), fixed = TRUE)

working_data <- t_lipid_data
working_data <- SID_to_metadata(working_data)
working_data <- character_to_factor(working_data) 


## summary biological & technical replicates
means_biol <- calc_by_replicate(working_data, "treatment", mean)
means_tech <- calc_by_replicate(working_data, "biol_replicate", mean)
means_all <- apply(dplyr::select_if(working_data, is.numeric), 2, mean)

sd_biol <- calc_by_replicate(working_data, "treatment", sd)
sd_tech <- calc_by_replicate(working_data, "biol_replicate", sd)
sd_all <- apply(dplyr::select_if(working_data, is.numeric), 2, sd)

var_biol <- calc_by_replicate(working_data, "treatment", var)
var_tech <- calc_by_replicate(working_data, "biol_replicate", var)
var_all <- apply(dplyr::select_if(working_data, is.numeric), 2, var)

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
## correlation heatmap --> for all variables
correlation_heatmap(working_data, interactive = TRUE)

### plots 
## paralell plot for <= 10 variables
parallel_plot(working_data, "treatment", plot_name) #todo Error in -groupColumn : invalid argument to unary operator 

## spider chart
spider_data <- SID_to_metadata(t_lipid_data) # calculate means so there is only one value per group
spider_data <- calc_by_replicate(spider_data, "treatment", mean)
rownames(spider_data) <- spider_data$Group.1

spider_chart(spider_data[1:10])


### PCA
## 1. check variances if scaling is necessary 
# (if there is a difference of > one potences between the variances)
wd <- working_data
groups <- wd$treatment
var_all <- apply(dplyr::select_if(wd, is.numeric), 2, var)

lipid_pca <- prcomp(select_if(wd[1:30], is.numeric), scale = TRUE, center = TRUE)
summary(lipid_pca)

var_lipid_pca <- lipid_pca$sdev ^ 2 # variance explained by each pc
prop_of_variance <- var_lipid_pca / sum(var_lipid_pca) # proportion of variance

scree_factoextra(lipid_pca)
scree_base(lipid_pca)

biplot_factoextra(lipid_pca, groups, ellipse = TRUE, loadings = FALSE)
biplot_ggplot2(input_df = wd, groups = "treatment", ellipse = TRUE, loadings = TRUE)

### Clustering

wd <- lipid_data
hc_wd <- hclust(dist(select_if(wd, is.numeric), "euclidean"), method = "average")

dir.create(paste(getwd(), "/examples", sep = ""), showWarnings = FALSE)
dir <- paste(getwd(), "/examples/lipComp", sep = "")

hclust_methods(wd)
hclust_performance_table(wd)
hclust_performance_plot(wd)
hclust_dendrogram(hc_wd)
hclust_heatmap(wd, row_names = wd$X)
hclust_heatmap_interactive(wd)


## hypothesis testing and volcano plot
wd_lps <- subset(working_data, working_data$treatment == "LPS")
wd_con <- subset(working_data, working_data$treatment == "Con")

t.test(wd_lps$`11-HETE_10.43`, wd_con$`11-HETE_10.43`)

