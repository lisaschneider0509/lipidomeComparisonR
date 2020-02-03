#### R script with functions for lipidome comparison ####

### Data handeling
## read and transpose data
read_transpose <- function(input_df){
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

## remove non numeric columns
# remove_non_numeric <- function(input_df){
  working_data[lapply(working_data, is.numeric)]
}

### Exploratory data analysis 
## calculate by biological replicate
# grouping by treatment calculates the mean/sd/etc. over all biological replicates
calc_biol_rep <- function(input_df, funct){ 
  as.data.frame(
    aggregate(dplyr::select_if(input_df, is.numeric), 
            by=list(input_df$treatment), 
            FUN=funct)
  )
}

## calculate by technical replicate
# grouping by biological replicate calculates the mean over all technical replicates
calc_tech_rep <- function(input_df, funct){
  as.data.frame(
    aggregate(dplyr::select_if(input_df, is.numeric), 
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
  input_df <- 
  pdf(paste(plot_name, "_boxplot", ".pdf", sep = ""))
  par(mfrow=c(3,3))
    for (i in 4:ncol(input_df[,1: ncol(input_df)])){
      col_name <- colnames(input_df)[i]
      boxplot(input_df[,i] ~ input_df[[by_factor]], 
              main = col_name,
              cex.main = 0.8,
              xlab = NULL, 
              ylab = NULL)}
  dev.off()}

## shapiro-wilk by factor
shapiro_by_factor <- function(input_df, by_factor){
  
  shapiro_statistic <- aggregate(dplyr::select_if(input_df, is.numeric), 
                                 by = list(input_df[[by_factor]]),
                                 FUN = function(x) {y <- shapiro.test(x); c(y$statistic)})
  
  shapiro_statistic <- tibble::add_column(shapiro_statistic, 
                                          value = "W", 
                                          .before = 1)
  
  shapiro_pvalue <- aggregate(dplyr::select_if(input_df, is.numeric), 
                              by = list(input_df[[by_factor]]),
                              FUN = function(x) {y <- shapiro.test(x); c(y$p.value)})
  
  shapiro_pvalue <- tibble::add_column(shapiro_pvalue, 
                                       value = "p-Value", 
                                       .before = 1)
  
  shapiro_all <- rbind(shapiro_statistic, shapiro_pvalue)
  shapiro_all[order(shapiro_all[,2]), ]
}

## test for correlation
correlation_plot <- function(input_df, method){
  panel.cor <- function(x, y){
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- round(cor(x, y, use = "pairwise.complete.obs", method = method), digits=2)
    txt <- paste0("R = ", r)
    text(0.5, 0.5, txt, cex = 1.5)
  }
  
  mypanel <- function(x, y){
    points(x, y, pch = 20)
  }
  
  pairs(dplyr::select_if(input_df, is.numeric),
        lower.panel = mypanel,
        upper.panel = panel.cor)
} # max 10 variables



