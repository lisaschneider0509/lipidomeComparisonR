### load packages
{
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
  library(ggpubr) # multiple plots on one page
  
  library(ggrepel)
  library(factoextra)
  library(ggfortify) # biplot with ggplot
  library(corrplot)
  library(FactoMineR)
  
  library(heatmaply) # interactive heatmap
  library(gplots) # heatmap
  library(plotly) # interactive ggplots
  library(htmlwidgets) # save plotly-plots as html
  library(dendextend)
  
  library(lipidr) # lipid set enrichment
  }

source("lipidome_comparison_dataTransformaions.R")
source("lipidome_comparison_EDA.R")
source("lipidome_comparison_pca.R")
source("lipidome_comparison_clustering.R")
source("lipidome_comparison_hypothesis_testing.R")

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

## set paths 
working_directory <- "/home/lisa/FH/Masterarbeit/LipidomeComparison"
setwd(working_directory)
meat_data_path <- "/home/lisa/FH/Masterarbeit/LipidomeComparison/data/meat_fish_final_raw.csv"
annotation_path <- "/home/lisa/FH/Masterarbeit/LipidomeComparison/data/meat_annotation.csv"

plot_path <- paste(working_directory, "/plots", sep = "")
plot_name <- paste(plot_path, "/meat_data", sep = "")

## meat data
### data processing
{
meat_data <- read.csv(meat_data_path, sep = ",", dec = ".", header = TRUE)
meat_data <- subset(meat_data, select = c(Compound, Type, Filename, Status, Group, Area))
meat_datameat_data <- subset(meat_data, Status == "Processed")
meat_data[meat_data=='N/F'] <- NA
meat_data$Area <- as.numeric(as.character(meat_data$Area))
meat_target <- subset(meat_data, Type == "Target Compound")
# meat_standard <- subset(meat_data, Type == "Internal Standard")

meat_target <- meat_data

meat_target <- subset(meat_target, select = c("Compound", "Filename", "Area"))

meat_target$Filename <- sub(".*probe","sample", meat_target$Filename)
meat_target$Filename <- sub("\\.*pos2","2", meat_target$Filename)
meat_target$Filename <- sub("\\.*pos","1", meat_target$Filename)

meat_target <- subset(meat_target, str_detect(meat_target$Filename, "sample") == TRUE)

meat_matrix <- flip_df_for_lipidr(meat_target)
meat_annotation <- read.csv(annotation_path)
}

tf_data <- read.csv(paste(getwd(), "/example_data/tracefinder_example.csv", sep = ""))
tf_data <- flip_df_for_lipidr(tf_data)
colnames(tf_data)[1] <- "lipids"
rownames(tf_data) <- tf_data$lipids
as_lipidomics_experiment(tf_data)


# annotation_data <- read.csv(paste(getwd(), "/example_data/tf_example_annotation2.csv", sep = ""))

# meat <- add_sample_annotation(meat, annotation_data)
plot_samples(meat, type = "boxplot", log = TRUE)
plot_lipidclass(meat, "boxplot")
meat_normalized <- normalize_pqn(meat, measure = "Area", exclude = "blank", log = TRUE)

mvaresults <- mva(meat_normalized, measure="Area", method="PCA")
plot_mva(mvaresults, color_by="group", components = c(1,2))
