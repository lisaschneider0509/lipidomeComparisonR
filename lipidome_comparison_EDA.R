### Exploratory data analysis  ###

#' Aggregate function results by factor
#' 
#' @description 
#' `calc_by_replicate` takes a data frame calculates the results of a fuction grouped by a factor
#' @details 
#' Take a data frame and calculate the results of a fuction grouped by a given factor. 
#' @param input_df a data frame with at least one factor column
#' @param factor a string with the column name to group by
#' @param funct a generic R function, that takes only one argument (e.g. mean(), summary(), etc.)
#' @examples 
#' calc_by_replicate(iris, "Species", mean)
calc_by_replicate <- function(input_df, factor, funct){ 
  as.data.frame(
    aggregate(dplyr::select_if(input_df, is.numeric), 
              by=list(input_df[[factor]]), 
              FUN=funct)
  )
}


### Graphical exploratory data analysis

#' Print one qqplot for each variable and group
#' 
#' @description `qqplot_by_factor` takes a data frame and prints a qq-plot for each group and varible
#' @details Take a data frame and prints a qq-plot for each group and varible. 
#' @param input_df a data frame with at least one factor column
#' @param factor a string with the column name to group by
#' @param out_path optional string. 
#' If out path is given, the generated plots are saved to a pdf document with the given name. 
#' @example 
#' qqplot_by_factor(iris, "Species")
#' \dontrun{
#' dir.create(paste(getwd(), "/examples", sep = ""), showWarnings = FALSE)
#' dir <- paste(getwd(), "/examples/iris", sep = "")
#' qqplot_by_factor(iris, "Species", dir)
#' }
qqplot_by_factor <- function(input_df, factor, out_path = "none"){
  levels <- levels(input_df[[factor]]) 
  numeric_df <- dplyr::select_if(input_df, is.numeric)
  
  func <- function(){
    par(mfrow=c(3,3))
    for (j in 1:ncol(numeric_df)){
      col_name <- colnames(numeric_df)[j]
      qqnorm(numeric_df[,j][input_df[[factor]] == levels[i]],
             main = paste(col_name, levels[i], sep = " "),
             cex.main = 0.8)
      qqline(numeric_df[,j][input_df[[factor]] == levels[i]])
    }
  }
  
  if(out_path == "none"){
    for (i in 1:length(levels)){ 
      func()
    }
  }
  else{
    
    for (i in 1:length(levels)){ 
      print(paste("Saving to ", out_path, "_qqplot_", levels[i], ".pdf", sep = ""))
      pdf(paste(out_path, "_qqplot_", levels[i], ".pdf", sep = ""))
      func()
      dev.off()
    }
  }
}

#' Print one histogram for each variable group. 
#' 
#' @description `histogram_by_factor` prints a histogram with density line
#' for each group and variable of a data frame
#' @details Generate a new grid/pdf dochument for each group of a data frame. 
#' Generate one plot for each variable of this subset data frame. 
#' @param input_df a data frame with at least one factor column
#' @param factor a string with the column name to group by
#' @param out_path optional string. 
#' @example 
#' histogram_by_factor(iris, "Species")
#' \dontrun{
#' dir.create(paste(getwd(), "/examples", sep = ""), showWarnings = FALSE)
#' dir <- paste(getwd(), "/examples/iris", sep = "")
#' histogram_by_factor(iris, "Species", dir)
#' }
histogram_by_factor <- function(input_df, factor, out_path = "none"){
  levels <- levels(input_df[[factor]]) 
  numeric_df <- dplyr::select_if(input_df, is.numeric)
  
  func <- function(){
    par(mfrow=c(3,3))
    for (j in 1:ncol(numeric_df)){
      col_name <- colnames(numeric_df)[j]
      hist(numeric_df[,j][input_df[[factor]] == levels[i]],
           main = paste(col_name, levels[i], sep = " "),
           cex.main = 0.8,
           xlab = NULL)
      lines(density(numeric_df[,j][input_df[[factor]] == levels[i]]))
      lines(density(numeric_df[,j][input_df[[factor]] == levels[i]], 
                    adjust = 1.5), lty = 2)
    }
  }
  
  if(out_path == "none"){
    for (i in 1:length(levels)){ 
      func()
    }
  }
  else{
    for (i in 1:length(levels)){ 
      print(paste("Saving to ", out_path, "_hist_", levels[i], ".pdf", sep = ""))
      pdf(paste(out_path, "_hist_", levels[i], ".pdf", sep = ""))
      func()
      dev.off()
    }
    
  }
}

#' Print boxplot by factor
#' 
#' @description `boxplot_by_factor` generates plot with one boxplot per group for each variable
#' @details A plot with one boxplot per group is generated for all variables of a data frame. 
#' @param input_df a data frame with at least one factor column
#' @param factor a string with the column name to group by
#' @param out_path optional string. 
#' @example 
#' boxplot_by_factor(iris, "Species")
#' \dontrun{
#' dir.create(paste(getwd(), "/examples", sep = ""), showWarnings = FALSE)
#' dir <- paste(getwd(), "/examples/iris", sep = "")
#' boxplot_by_factor(iris, "Species", dir)
#' }
#'
boxplot_by_factor <- function(input_df, factor, out_path = "none"){
  numeric_df <- dplyr::select_if(input_df, is.numeric)
  
  func <- function(){
    for (i in 1:ncol(numeric_df)){
      col_name <- colnames(input_df)[i]
      boxplot(numeric_df[,i] ~ input_df[[factor]],
              main = col_name,
              xlab = NULL,
              ylab = NULL)
    }
  }
  
  if(out_path == "none"){
    par(mfrow=c(3,3), 
        cex.main = 1, 
        cex.axis = 0.7)
    func()
  }
  else{
    print(paste("Saving to ", out_path, "_boxplot.pdf", sep = ""))
    pdf(paste(out_path, "_boxplot.pdf", sep = ""))
    par(mfrow=c(3,3), 
        cex.main = 1, 
        cex.axis = 0.7)
    func()
    dev.off()
  }
}


#' Assess normality for each group
#' 
#' @description `shapiro_by_factor` takes a data frame and applies the shapiro-wilk test 
#' to every group and variable
#' @details Shapiro wilk test is applied to every group and variable of a data frame. 
#' The results are aggregated and printed in a table. 
#' @param input_df a data frame with at least one factor column
#' @param factor a string with the column name to group by
#' @param out_path optional string. 
#' @example 
#' shapiro_by_factor(iris, "Species")
#' \dontrun{
#' dir.create(paste(getwd(), "/examples", sep = ""), showWarnings = FALSE)
#' dir <- paste(getwd(), "/examples/iris", sep = "")
#' shapiro_by_factor(iris, "Species", dir)
#' }
shapiro_by_factor <- function(input_df, factor, out_path = "none"){
  
  
  shapiro_statistic <- aggregate(dplyr::select_if(input_df, is.numeric), 
                                 by = list(input_df[[factor]]),
                                 FUN = function(x) {y <- shapiro.test(x); c(y$statistic)})
  
  shapiro_statistic <- tibble::add_column(shapiro_statistic, 
                                          value = "W", 
                                          .before = 1)
  
  shapiro_pvalue <- aggregate(dplyr::select_if(input_df, is.numeric), 
                              by = list(input_df[[factor]]),
                              FUN = function(x) {y <- shapiro.test(x); c(y$p.value)})
  
  shapiro_pvalue <- tibble::add_column(shapiro_pvalue, 
                                       value = "p-Value", 
                                       .before = 1)
  
  shapiro_all <- rbind(shapiro_statistic, shapiro_pvalue)
  shapiro_all <- shapiro_all[order(shapiro_all[,2]), ]
  
  if(out_path == "none"){
    DT::datatable(shapiro_all)
  } else {
    DT::datatable(shapiro_all)
    print(paste("Writing to", out_path, "_shapiro.csv", sep = ""))
    write.csv(shapiro_all, paste(out_path, "_shapiro.csv", sep = ""))
  }
}


#' Correlation plots 
#' 
#' @description `correlation_plot` calculates correlations of variables and displays them in a dotplot
#' @details A grid of dotplots displaying the correlations between the variables of a data frame is generated. 
#' Use with less than 10 variables ideally. 
#' @param input_df data frame. 
#' @param method string. Method for calculating the correlation. 
#' Options: "pearson", "kendall", "spearman" (default). 
#' @param out_path optional string. Path to save correlation plot to png. 
#' If out_path is empty the correlation plot is printed to the device. 
#' @example 
#' correlation_plot(iris)
#' \dontrun{
#' dir.create(paste(getwd(), "/examples", sep = ""), showWarnings = FALSE)
#' dir <- paste(getwd(), "/examples/iris", sep = "")
#' correlation_plot(iris, out_path = dir)
#' }
correlation_plot <- function(input_df, method = "spearman", out_path = "none"){
  out_name <- paste(out_path, "_correlations", ".png", sep = "")
  
  func <- function(){
    par(cex.main = 1, cex.axis = 0.8, 
        family = "sans", font = 1)
    colour = viridis(n = 2)
    pairs.panels(select_if(input_df, is.numeric),
                 method = method,
                 hist.col = colour[2], rug = FALSE,
                 density = TRUE,  lm = TRUE, ci = FALSE, col = colour[1],
                 ellipses = FALSE,
                 pch = 20, cex = 0.7,
                 cex.cor = 0.5, main = "Correlation plot")
  }
  
  if(out_path != "none"){
    print(paste("Saving to ", out_name, sep = ""))
    png(out_name)
    func()
    dev.off() 
    func()
  }
  else{
    func()
  }
}

#' Correlation heatmap 
#' 
#' @description `correlation_heatmap` calculates correlations of variables and displays them in a heatmap
#' @details A heatmap displaying the correlations between the variables of a data frame is generated. 
#' The heatmap is optionally interactive. 
#' @param input_df data frame. 
#' @param method string. Method for calculating the correlation. 
#' Options: "pearson", "kendall", "spearman" (default). 
#' @param interactive logical. Print heatmap to device (FALSE, default) 
#' or open interactive heatmap in browser (TRUE). 
#' @param out_path string. Path to save heatmap to png. 
#' If out_path is empty the heatmap is printed to the device. 
#' @example 
#' correlation_heatmap(iris)
#' correlation_heatmap(iris, interactive = TRUE)
#' \dontrun{
#' dir.create(paste(getwd(), "/examples", sep = ""), showWarnings = FALSE)
#' dir <- paste(getwd(), "/examples/iris", sep = "")
#' correlation_heatmap(iris, interactive = TRUE, out_path = dir)
#' correlation_heatmap(iris, out_path = dir)
#' }
correlation_heatmap <- function(input_df, 
                                method = "spearman", 
                                interactive = FALSE, 
                                out_path = "none"){
  cor_matrix <- cor(dplyr::select_if(input_df, is.numeric), method = method)
  melted_cor_matrix <- reshape::melt(cor_matrix)
  names(melted_cor_matrix) <- c("x", "y", "correlation")
  head(melted_cor_matrix)
  
  cor_heatmap <- ggplot(data = melted_cor_matrix, aes(x=x, y=y, fill=correlation)) +
    geom_tile() +
    ggtitle("Spearman correlation") +
    scale_fill_viridis_c(option = "magma") +
    my_theme +
    theme(axis.title = element_blank(), 
          axis.text.y = element_text(size = 6),
          axis.text.x = element_text(angle = 90, size = 6, hjust = 1))
  
  if(out_path != "none"){
    print(paste("Saving heatmap to ", out_path, "_cor_heatmap.png", sep = ""))
    ggsave(paste(out_path, "_cor_heatmap.png", sep = ""), 
           plot = cor_heatmap)
  }
  
  if (interactive == TRUE) {
    plotly::ggplotly(cor_heatmap) # interactive heatmap
  } else {
    cor_heatmap # static heatmap
  }
  
}

#' Parallel coordinates plot
#' 
#' @description `parallel_plot` prints a paralell coordinates plot using a data frame
#' @details This function takes a data frame with at least one factor variable and 
#' displays #' it in a pralell coordinates plot, where the different groups are 
#' color coded. This plot works best for <= 10 parameters. 
#' @param input_df data frame. 
#' @param groupColumn numeric. Column to sort by.  
#' @param out_path string. Path to save parallel plot to png. 
#' If out_path is empty the parallel plot is printed to the device.
#' @param titles a vector of strings. 
#' 1. Main Title. Default = "Parallel plot"
#' 2. X-axis title. Default = ""
#' 3. Y-axis title. Default = ""
#' @param scale string. Method used to scale the variable. Default = "globalminmax". 
#' Options: "std", "robust", "uniminmax", "globalminmax", "center", "centerObs". 
#' For more information on the options see help(ggparcoord). 
#' @example 
#' parallel_plot(iris, 5)
#' parallel_plot(iris, 5, scale = "center", titles = c("Centered parallel plot", "Parameters", "Univariate scale to standardize vertical height"))
#' \dontrun
#' dir.create(paste(getwd(), "/examples", sep = ""), showWarnings = FALSE)
#' dir <- paste(getwd(), "/examples/iris", sep = "")
#' parallel_plot(iris, 5, out_path = dir)
parallel_plot <- function(input_df,  groupColumn, out_path = "none", 
                          titles = c("Parallel Plot", "", ""), 
                          scale = "globalminmax"){
  parallel_title <- titles[[1]]
  x_axis_title <- titles[[2]]
  y_axis_title <- titles[[3]]
  x_labels <- colnames(input_df[-groupColumn])
  y_labels <- NULL
  cols <- c(1:ncol(input_df))
  
  par_plot <- ggparcoord(input_df,
                         columns = cols[-groupColumn],
                         groupColumn = groupColumn,
                         showPoints = TRUE,
                         scale = scale, # "globalminmax" is default
                         alphaLines = 0.5
  )  +
    ggtitle(parallel_title) +
    xlab(x_axis_title) +
    ylab(y_axis_title) +
    scale_color_viridis(discrete=TRUE) +
    scale_x_discrete(breaks = colnames(input_df[-groupColumn]),
                     labels = x_labels) +
    geom_point(shape = 20, size = 0.5) +
    my_theme
  
  if(out_path != "none"){
    print(paste("Saving parallel coordinates plot to ", out_path, "_parcoord.png", sep = ""))
    ggsave(paste(out_path, "_parcoord.png", sep = ""),
           plot = par_plot)
  }
  else{
    par_plot
  }
}


#' Spider Chart
#' 
#' @description `spider chart` takes a minimized data frame (one value per group) and prints a spider chart
#' @details This function takes a data frame of with one value per group (i.e. calculate mean groupwise). 
#' This minimized data frame is used to draw a spider chart (also radar chart odr network plot). 
#' The ideal numer of parameters for a spider chart is <= 10. Also the shape of the graph depends on 
#' the order of parameters. If the data has large differences in size, normalizing or scaling the data is necessary.
#' Doesn't work with non normal data.  
#' @param minimized_df data frame. A data frame that contains only one value per group and variable. 
#' Most often this will be a data frame of means calculated from another data frame. 
#' @param tile string. Main title of the chart. Default = "Spider chart"
#' @param out_path string. Path to save spider chart to png. 
#' If out_path is empty, the spider chart is printed to the device.
#' @example 
#' minimized_iris <- aggregate(dplyr::select_if(iris, is.numeric), 
#'                             by = list(iris$Species), 
#'                             FUN = mean)
#' rownames(minimized_iris) <- minimized_iris$Group.1
#' spider_chart(minimized_iris, title = "Spider chart of iris")
#' \dontrun
#' minimized_iris <- aggregate(dplyr::select_if(iris, is.numeric), 
#'                             by = list(iris$Species), 
#'                             FUN = mean)
#' rownames(minimized_iris) <- minimized_iris$Group.1
#' dir.create(paste(getwd(), "/examples", sep = ""), showWarnings = FALSE)
#' dir <- paste(getwd(), "/examples/iris", sep = "")
#' spider_chart(minimized_iris, out_path = dir)
spider_chart <- function(minimized_df, title="Spider chart", out_path = "none"){ # todo get labels ot of the plot
  
  out_name <- paste(out_path, "_spiderChart", ".png", sep = "")
  
  func <- function(){
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
               centerzero = FALSE)
    
    title(main = title, cex.main = 0.9, font.main = 1)
    
    ## Add a legend
    legend(x=-2, 
           y=1.1, 
           legend = rownames(spider_data[-(1:2),]), 
           bty = "n", 
           pch=20, 
           col=colors_border, 
           text.col = "black", 
           cex=0.7, pt.cex=1.3)
  }
  
  ## print to device or save
  if(out_path != "none"){
    print(paste("Saving to ", out_name))
    png(filename = out_name)
    func()
    dev.off()
  }
  else{
    func()
  }
}


