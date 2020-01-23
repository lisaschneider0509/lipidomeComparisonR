### Install packages
# install.packages("ggplot2")
# install.packages("stringr")
install.packages("tidyr")
install.packages("data.table")
install.packages("textshape")
install.packages("tibble")

### load packages
library(gridExtra)
library(stringr)
library(ggplot2)
library(tidyr)
library(data.table)
library(textshape)
library(tibble)

theme_set(
  theme_minimal() +
    theme(legend.position = "top")
)

# set variables
working_directory <- "/home/lisa/FH/Masterarbeit/LipidomeComparison"
setwd(working_directory)

input_path <- "/home/lisa/FH/Masterarbeit/LipidomeComparison/data/Probe-Datensatz_lisa.csv"
test_path <- "/home/lisa/FH/Masterarbeit/LipidomeComparison/data/prefix_data.csv"
plot_path <- paste(working_directory, "/plots", sep = "")
plot_name <- paste(plot_path, "/test_data", sep = "")


if (! file.exists(plot_path)){
  print("HI")
} else {
  print("Nope")
}

# load & transform data
lipid_data <- read.csv(input_path, sep = ",", dec = ".", header = TRUE) #read data
t_lipid_data <- transpose(as.data.table(lipid_data), make.names = 1)
t_lipid_data <- as.data.frame(t_lipid_data)
row.names(t_lipid_data) <- colnames(lipid_data[-1])


test_data <- read.csv(test_path, sep = ",", dec = ".", header = TRUE) #read data
t_test_data <- transpose(as.data.table(test_data), make.names = 1)
t_test_data <- as.data.frame(t_test_data)
row.names(t_test_data) <- colnames(test_data[-1])

working_data <- t_test_data



### summary biological & technical replicates
meta_info <- read.table(text = row.names(working_data), sep = "_")
meta_info$V2 <- paste(meta_info$V1, meta_info$V2, sep = "_")
meta_info$V3 <- NULL
working_data <- tibble::add_column(working_data, treatment = meta_info$V1, .before = 1)
working_data <- tibble::add_column(working_data, biol_replicate = meta_info$V2, .after = 1)

means_biol <- as.data.frame(aggregate(working_data[-(1:2)], by=list(working_data$treatment), FUN=mean)) # grouping by sample-ID calculates the mean over all biological replicates
means_tech <- as.data.frame(aggregate(working_data[-(1:2)], by=list(working_data$biol_replicate), FUN=mean)) # grouping by biological replicate calculates the mean over all technical replicates

sd_biol <- as.data.frame(aggregate(working_data[-(1:2)], by=list(working_data$treatment), FUN=sd)) # grouping by sample-ID calculates the mean over all biological replicates
sd_tech <- as.data.frame(aggregate(working_data[-(1:2)], by=list(working_data$biol_replicate), FUN=sd)) # grouping by biological replicate calculates the mean over all technical replicates

## plots for normal distribution
sample_df <- working_data

pdf(paste(plot_name, "_qqplot", ".pdf", sep = ""))
par(mfrow=c(3,3))
for (i in 3:ncol(sample_df[,1: ncol(sample_df)])){
  col_name <- colnames(sample_df)[i]
  levels <- levels(sample_df$treatment)
  color_list <- c("blue", "red", "green") 
  
  # qqnorm(sample_df[,i][sample_df$treatment == levels[1]], main = col_name, cex.main = 0.8, xlab = NULL, col = color_list[1])
  # qqline(sample_df[,i][sample_df$treatment == levels[1]])
  # 
  # for (j in 2:length(levels)){
  #   qqnorm(sample_df[,i][sample_df$treatment == levels[j]], main = col_name, cex.main = 0.8, xlab = NULL, col = color_list[j], add = TRUE)
  #   qqline(sample_df[,i][sample_df$treatment == levels[j]], col = color_list[j], add = TRUE)
  # }
  
  qqnorm(sample_df[,i], main = col_name, cex.main = 0.8) # qq-plot
  qqline(sample_df[,i]) # expected values from qq if data was ND
}
dev.off()
 
pdf(paste(plot_name, "_histogram", ".pdf", sep = ""))
par(mfrow=c(3,3))
for (i in 3:ncol(sample_df[,1: ncol(sample_df)])){
  col_name <- colnames(sample_df)[i]
  levels <- levels(sample_df$treatment)
  color_list <- c(rgb(173,216,230,max = 255, alpha = 80, names = "lt.blue"), rgb(255,192,203, max = 255, alpha = 80, names = "lt.pink"), rgb(144,238,144, max = 255, alpha = 80, names = "lt.green")) 

  # hist(sample_df[,i][sample_df$treatment == levels[1]], main = col_name, cex.main = 0.8, xlab = NULL, col = color_list[1], border = FALSE)
  # lines(density(sample_df[,i][sample_df$treatment == levels[1]]))
  # 
  # for (j in 2:length(levels)){
  #   hist(sample_df[,i][sample_df$treatment == levels[j]], add = TRUE, col = color_list[j], border = FALSE)
  # }
  hist(sample_df[,i], main = col_name, xlab = "intensity", cex.main = 0.8)
  lines(density(sample_df[,i]))
  lines(density(sample_df[,i],adjust=1.5),col=2) # 1.5 x bandwidth
}
dev.off()

## test for normal distribution
### Don't use with multi modal data --> check histogram and qq plots first




