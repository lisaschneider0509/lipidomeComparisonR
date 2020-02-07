#### R script with functions for lipidome comparison ####

### Data handeling
## read and transpose data
pretty_transpose <- function(input_df){
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

## save all non-numeric columns of a dataframe as factor
character_to_factor <- function(input_df){
  for(i in 1:ncol(input_df)){
    if (!is.numeric(input_df[,i])){
      input_df[,i] <- factor(input_df[,i])
    }
  }
  new_df <- input_df
}

## remove non numeric columns
# remove_non_numeric <- function(input_df){
#   working_data[lapply(working_data, is.numeric)]
# }

### Exploratory data analysis 
## calculate by biological replicate
# grouping by treatment calculates the mean/sd/etc. over all biological replicates

calc_by_replicate <- function(input_df, by_factor, funct){ 
  as.data.frame(
    aggregate(dplyr::select_if(input_df, is.numeric), 
              by=list(input_df[[by_factor]]), 
              FUN=funct)
  )
}

### Graphical exploratory data analysis 
## qqplot by factor
# print one pdf document with qq-plotsper factor for all lipids
qqplot_by_factor <- function(input_df, by_factor, out_path){
  levels <- levels(input_df[[by_factor]])
  for (i in 1:length(levels)){
    # pdf(paste(out_path, "_qqplot_", levels[i], ".pdf", sep = ""))
    par(mfrow=c(3,3))
    for (j in 4:ncol(input_df[,1: ncol(input_df)])){
      col_name <- colnames(input_df)[j]
      
      qqnorm(input_df[,j][input_df[[by_factor]] == levels[i]],
             main = paste(col_name, levels[i], sep = " "),
             cex.main = 0.8)
      qqline(input_df[,j][input_df[[by_factor]] == levels[i]])
    }
    # dev.off()
  }
}

## histogram by factor
# print one pdf document with histograms and density lines per factor for all lipids
histogram_by_factor <- function(input_df, by_factor, out_path){
  levels <- levels(input_df[[by_factor]])
  for (i in 1:length(levels)){
    # pdf(paste(plot_name, "_histogram_", levels[i], ".pdf", sep = ""))
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
    # dev.off()
  }
}

## boxplot by factor
boxplot_by_factor <- function(input_df, by_factor, out_path){ 
  # pdf(paste(plot_name, "_boxplot", ".pdf", sep = ""))
  par(mfrow=c(3,3))
    for (i in 4:ncol(input_df[,1: ncol(input_df)])){
      col_name <- colnames(input_df)[i]
      boxplot(input_df[,i] ~ input_df[[by_factor]], 
              main = col_name,
              cex.main = 0.8,
              xlab = NULL, 
              ylab = NULL)}
  # dev.off()
}

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
  # pdf(paste(plot_name, "_correlations", ".pdf", sep = ""))
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
  # dev.off()
} # max 10 variables

## parallel coordinates plot 
parallel_plot <- function(input_df,  factor, out_path, 
                          titles = c("Paralell Plot", "", ""), 
                          scale = "center"){
  # # pdf(paste(plot_name, "_ParallelPlot", ".pdf", sep = ""))
  par(mfrow=c(1,1))
  
  parallel_title <- titles[[1]]
  x_axis_title <- titles[[2]]
  y_axis_title <- titles[[3]]
  x_labels <- substring(colnames(working_data),
                        first = 1,
                        last = 8)
  y_labels <- NULL
  
  # Plot
  ggparcoord(data,
             columns = 1:10, 
             groupColumn = ncol(data),
             showPoints = TRUE, 
             scale="center", # "center" is default
             alphaLines = 0.3) + 
    ggtitle(parallel_title) +
    xlab(x_axis_title) + 
    ylab(y_axis_title) +
    scale_color_viridis(discrete=TRUE) +
    theme_ipsum()+
    theme( # font settings
      plot.title = element_text(size=14, hjust = 0.5), 
      axis.text.x = element_text(size = 8), 
      axis.title = element_text(size = 10), 
      axis.title.x = element_text(size = 10, hjust = 0.5),
      axis.title.y = element_text(size = 10, hjust = 0.5),
      legend.text = element_text(size = 8), 
      legend.title = element_text(size = 10)) + 
    scale_x_discrete(breaks = colnames(working_data), 
                     labels = my_labels)
  
  # # dev.off()
}

## spider chart (= radar chart, network plot, etc.)
spider_chart <- function(minimized_df, title="Spider chart", out_path){ # todo get labels ot of the plot
  ## input_df <= 10 columns 
  ## minimized_df = dataframe with only one row per group (i.e. calculate means)
  
  # pdf(paste(plot_name, "_ParallelPlot", ".pdf", sep = ""))
  spider_legend <- row.names(minimized_df) # set new row names 
  spider_data <- dplyr::select_if(minimized_df, is.numeric) # remove column with rownames
  spider_labels <- substring(colnames(spider_data), first = 1, last = 6) # set max. label length to 10 characters

  spider_min <- floor(min(spider_data))
  spider_max <- ceiling(max(spider_data))
  spider_data <- as.data.frame(select_if(spider_data, is.numeric))

  # add max and min to the dataframe to plot the grid
  spider_data <- rbind(spider_min, spider_max, spider_data)

  ## set colors
  colors_border = as.vector(viridis(n = nrow(minimized_df), option = "viridis"))
  colors_in = alpha(colors_border, alpha = 0.1)

  ## radar chart
  radarchart(spider_data,
              axistype=0,
              #custom polygon
              pcol=colors_border,
              pfcol=colors_in,
              plwd=4,
              plty=1,
              # custom grid
              cglcol="grey",
              cglty=1,
              axislabcol="grey",
              cglwd=0.8,
              # custom labels
              vlcex=0.6, 
             centerzero = FALSE, 
             title = title
  )

  ## Add a legend
  legend(x=-2, 
         y=1.1, 
         legend = rownames(spider_data[-(1:2),]), 
         bty = "n", 
         pch=20, 
         col=colors_border, 
         text.col = "black", 
         cex=0.7, pt.cex=1.3)
  # dev.off()
}

