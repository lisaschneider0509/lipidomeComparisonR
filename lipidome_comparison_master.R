### load packages
library(dplyr) # select part of data
library(stringr) # count separators
library(data.table) # transpose data frame
library(impute)
library(imputeLCMD)

library(ggplot2)#, # plots
library(tibble) # data frame manipulation
library(viridis) # colorblind save color schemes
library(GGally) # paralell plot
library(fmsb) # spider chart
library(scales) # scale opacity of filling (alpha)

library(ggrepel)
library(factoextra)
library(ggfortify) # biplot with ggplot

library(heatmaply) # interactive heatmap
library(gplots) # heatmap
library(plotly) # interactive heatmap

# library(psych) # for correlation plot 
# library(gridExtra)
# library(devtools)
# library(dendextend)
# library(limma) # hypothesis testing


source("lipidome_comparison_dataTransformaions.R")
source("lipidome_comparison_EDA.R")
source("lipidome_comparison_pca.R")
source("lipidome_comparison_clustering.R")
source("lipidome_comparison_hypothesis_testing.R")

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

test_path <- "/home/lisa/FH/Masterarbeit/LipidomeComparison/data/Probe-Datensatz_lisa.csv"
meat_data_path <- "/home/lisa/FH/Masterarbeit/LipidomeComparison/data/meat_fish_final_raw.csv"

plot_path <- paste(working_directory, "/plots", sep = "")
plot_name <- paste(plot_path, "/meat_data", sep = "")

## meat data
### data processing

meat_data <- read.csv(meat_data_path, sep = ",", dec = ".", header = TRUE)
meat_data <- subset(meat_data, select = c(Compound, Type, Filename, Status, Group, Area))
meat_data <- subset(meat_data, Status == "Processed")
meat_data[meat_data==''] <- NA
meat_data[meat_data=='N/F'] <- NA
meat_data$Area <- as.numeric(meat_data$Area)
meat_target <- subset(meat_data, Type == "Target Compound")
meat_standard <- subset(meat_data, Type == "Internal Standard")

meat_target <- flip_df(meat_target)


meat_target <- subset(meat_target, !is.na(Group))
meat_target$SID <- sub(".*probe","sample", meat_target$SID)
meat_target$SID <- sub("\\.*pos2","2", meat_target$SID)
meat_target$SID <- sub("\\.*pos","1", meat_target$SID)

meat_AS <- meat_target$SID[str_detect(meat_target$SID, "AS") == TRUE] 
meat_target$SID[str_detect(meat_target$SID, "AS") == TRUE] <- sub(".*sample","AS_sample", meat_AS)
meat_N <- meat_target$SID[str_detect(meat_target$SID, "AS") == FALSE] 
meat_target$SID[str_detect(meat_target$SID, "AS") == FALSE] <- sub(".*sample","N_sample", meat_N)
meat_target$SID <- str_remove(meat_target$SID, "_AS")

meta_info <- read.table(text = meat_target$SID, sep = "_")
colnames(meta_info) <- c("Treatment", "Sample_nr", "Biol_rep", "Tech_rep")
meta_info$Biol_rep <- paste(meta_info$Sample_nr, meta_info$Biol_rep, sep = "_")
meta_info$Tech_rep <- paste(meta_info$Biol_rep, meta_info$Tech_rep, sep = "_")

meat_target <- cbind(meat_target$SID, meta_info, meat_target[, -1])
meat_target <- droplevels(meat_target)
levels(meat_target$Group)[levels(meat_target$Group) == "fleisch"] <- "meat"
levels(meat_target$Group)[levels(meat_target$Group) == "wild"] <- "game"
levels(meat_target$Group)[levels(meat_target$Group) == "FISCH"] <- "fish"
colnames(meat_target) <- c("SID", colnames(meat_target[-1]))
meat_N <- subset(meat_target, Treatment == "N")
meat_AS <- subset(meat_target, Treatment == "AS")


levels(meat_N$Group)

## Exploratory data analysis

### impute missing values #todo find option to avoid imputation of negative values
#### remove columns where all values are missing
impute_meat <- meat_N[, which(colMeans(!is.na(meat_N)) > 0.8)] 
impute_meat <- as.matrix(select_if(impute_meat, is.numeric))

#### perform missing data imputation
meat_QRILC <- impute.QRILC(impute_meat, tune.sigma = 1)
# meat_MinDet <- impute.MinDet(impute_meat)

# meat_imputed <- meat_MinDet
meat_imputed <- as.data.frame(meat_QRILC[[1]])
meat_imputed <- cbind(meat_N[, 1:6], meat_imputed)
meat_imputed <- droplevels(meat_imputed) # remove unused levels from factors

#### calculate the means for the replicates
{
meat_groups <- generate_categorical_table(meat_imputed$Group)
meat_treatment <- generate_categorical_table(meat_imputed$Treatment)

meat_numeric <- meat_imputed
meat_numeric$Group <- as.numeric(meat_numeric$Group)
meat_biol <- calc_by_replicate(meat_numeric, meat_numeric$Sample_nr, mean)
meat_tech <- calc_by_replicate(meat_numeric, meat_numeric$Biol_rep, mean)

nmb <- paste_catecorical_variable(meat_biol, 2, meat_groups)
nmt <- paste_catecorical_variable(meat_tech, 2, meat_groups)
}

### graphical exploratory data analysis
# qqplot_by_factor(meat_imputed, "Group", out_path = plot_name)
# histogram_by_factor(meat_imputed, "Group", out_path = plot_name)
# boxplot_by_factor(meat_imputed, "Group", out_path = plot_name)
# 
# parallel_plot(meat_imputed, meat_imputed$Group, out_path = plot_name)
# meat_spider <- calc_by_replicate(meat_imputed, meat_imputed$Group, mean)
# spider_chart(meat_spider, legend_lab = meat_spider$Group.1, out_path = plot_name)

### test for normality
meat_normality <- shapiro_by_factor(meat_imputed, meat_imputed$Group)

### test for correlation
meat_correlation <- cor(select_if(meat_imputed, is.numeric), method = "spearman")
correlation_heatmap(meat_imputed, interactive = TRUE, out_path = plot_name)

### PCA
meat_pca <- prcomp(select_if(meat_imputed, is.numeric))

meat_pca_var <- meat_pca$sdev ^ 2
prop_var_meat <- round(meat_pca_var / sum(meat_pca_var) * 100, 2)
cum_prop_var_meat <- cumsum(prop_var_meat * 100)
proportion_of_variance_table <- data.frame(Proportion_of_variance = prop_var_meat, Cummulative_proportion_of_variance = cum_prop_var_meat)

scree_base(meat_pca)
scree_factoextra(meat_pca)
biplot_ggplot2(meat_imputed, "Group", loadings = FALSE, ellipse = TRUE)
biplot_factoextra(meat_pca, meat_imputed$Group, ellipse = TRUE)

### Clustering 
meat_clust <- data.frame(Group = meat_imputed$Group)
meat_clust <- cbind(meat_clust, select_if(meat_imputed, is.numeric))
rownames(meat_clust) <- meat_imputed$SID

hclust_performance_table(meat_clust)
hclust_performance_plot(meat_clust)

meat_dist <- dist(select_if(meat_clust, is.numeric), method = "manhattan")
meat_hclust <- hclust(meat_dist, method = "average")
hclust_dendrogram(meat_hclust, labs = meat_clust$Group)

hclust_heatmap(meat_clust, dist_method = "manhattan", hclust_method = "average", row_names = meat_clust$Group)
hclust_heatmap_interactive(meat_clust, dist_method = "manhattan", hclust_method = "average")

### hypothesis testing & volcano plot
meat_vs_fish <- subset(meat_imputed, Group == "fish" | Group == "meat")
meat_vs_fish <- droplevels(meat_vs_fish)

p_meat_vs_fish <- one_sample_test_by_col(meat_vs_fish, meat_vs_fish$Group, method = wilcox.test)
adj_meat_vs_fish <- p.adjust(p_meat_vs_fish$p_values, method = "fdr")
fc_meat_vs_fish <- log2_foldchange(meat_vs_fish, meat_vs_fish$Group)

meat_volcano <- data.frame(p_value = p_meat_vs_fish, adj_p_value = adj_meat_vs_fish, log2_foldchange = fc_meat_vs_fish)
# meat_volcano$log2_foldchange[is.nan(meat_volcano$log2_foldchange)] <- NA
meat_volcano <- meat_volcano[complete.cases(meat_volcano),]
# meat_volcano$log2_foldchange[is.na(meat_volcano$log2_foldchange)] <- median(meat_volcano$log2_foldchange, na.rm = TRUE)


volcano_plot(meat_volcano, 
             foldchange_col = meat_volcano$log2_foldchange, 
             significance_col = meat_volcano$p_values, 
             foldchange = 2, 
             significance = 0.5)


