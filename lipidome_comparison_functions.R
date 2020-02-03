#### R script with functions for lipidome comparison ####

### Data handeling
## read and transpose data
read_transpose <- function(input_path){
  input_df <- read.csv(input_path, sep = ",", dec = ".", header = TRUE) #read data
  t_input_df <- data.table::transpose(as.data.table(input_df), make.names = 1)
  t_input_df <- as.data.frame(t_input_df)
  row.names(t_input_df) <- colnames(input_df[-1])
  t_input_df
}

## get meta information from sample name
SID_to_metadata <- function(input_df){
  n_sep <- mean(stringr::str_count(row.names(input_df), "_"))
  meta_info <- read.table(text = row.names(input_df), sep = "_")
  
  if(n_sep == 2){
    biol_rep <- paste(meta_info$V1, meta_info$V2, sep = "_")
    treatment <- meta_info$V1
    input_df <- tibble::add_column(input_df, biol_replicate = biol_rep, .before = 1)
    input_df <- tibble::add_column(input_df, treatment = treatment, .before = 1)
    input_df <- tibble::add_column(input_df, experiment = NA, .before = 1)
  } else if(n_sep == 3){
    biol_rep <- paste(meta_info$V1, meta_info$V2, meta_info$V3, sep = "_")
    treatment <- meta_info$V2
    experiment <- meta_info$V1
    input_df <- tibble::add_column(input_df, biol_replicate = biol_rep, .before = 1)
    input_df <- tibble::add_column(input_df, treatment = treatment, .before = 1)
    input_df <- tibble::add_column(input_df, experiment = experiment, .before = 1)
  } else{
    print("Wrong format")
  }
}

### Exploratory data analysis 
## calculate by biological replicate
# grouping by treatment calculates the mean/sd/etc. over all biological replicates
calc_biol_rep <- function(input_df, funct){ 
  as.data.frame(
    aggregate(input_df[-(1:3)], 
            by=list(input_df$treatment), 
            FUN=funct)
  )
}

## calculate by technical replicate
# grouping by biological replicate calculates the mean over all technical replicates
calc_tech_rep <- function(input_df, funct){
  as.data.frame(
    aggregate(input_df[-(1:3)], 
              by=list(input_df$biol_replicate), 
              FUN=funct)
  )
}
  
### Graphical exploratory data analysis 
## qqplot by factor
# print one pdf document with qq-plotsper factor for all lipids
qqplot_by_factor <- function(input_df, by_factor, out_path){
  levels <- levels(input_df[[by_factor]])
  for (i in 1:length(levels)){
    pdf(paste(out_path, "_qqplot_", levels[i], ".pdf", sep = ""))
    par(mfrow=c(3,3))
    for (j in 4:ncol(input_df[,1: ncol(input_df)])){
      col_name <- colnames(input_df)[j]
      
      qqnorm(input_df[,j][input_df[[by_factor]] == levels[i]],
             main = paste(col_name, levels[i], sep = " "),
             cex.main = 0.8)
      qqline(input_df[,j][input_df[[by_factor]] == levels[i]])
    }
    dev.off()
  }
}

## histogram by factor
# print one pdf document with histograms and density lines per factor for all lipids
histogram_by_factor <- function(input_df, by_factor, out_path){
  levels <- levels(input_df[[by_factor]])
  for (i in 1:length(levels)){
    pdf(paste(plot_name, "_histogram_", levels[i], ".pdf", sep = ""))
    par(mfrow=c(3,3))
    for (j in 4:ncol(input_df[,1: ncol(input_df)])){
      col_name <- colnames(input_df)[j]
      
      hist(input_df[,j][input_df[[by_factor]] == levels[i]],
           main = col_name,
           cex.main = 0.8,
           xlab = NULL)
      lines(density(input_df[,j][input_df[[by_factor]] == levels[i]]))
      lines(density(input_df[,j][input_df[[by_factor]] == levels[i]], adjust = 1.5), lty = 2)
    }
    dev.off()
  }
}

## boxplot by factor
boxplot_by_factor <- function(input_df, by_factor, out_path){
    pdf(paste(plot_name, "_boxplot", ".pdf", sep = ""))
  par(mfrow=c(3,3))
    for (i in 4:ncol(input_df[,1: ncol(input_df)])){
      col_name <- colnames(input_df)[i]
      boxplot(working_data[,i] ~ working_data[[by_factor]], 
              main = col_name,
              cex.main = 0.8,
              xlab = NULL, 
              ylab = NULL)}
  dev.off()}

