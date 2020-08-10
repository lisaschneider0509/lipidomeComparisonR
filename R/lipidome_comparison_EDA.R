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
#' @param na_action function. Indcates what should happen if there are NA values in the data. c(NULL (default), na.omit).
#' @examples 
#' calc_by_replicate(iris, iris$Species, mean)
calc_by_replicate <- function(input_df, 
                              factor, 
                              funct, 
                              na_action = na.omit
                              ){ 
  as.data.frame(
    aggregate(select_if(input_df, is.numeric), 
              by=list(factor), 
              FUN=funct, 
              na.action = na_action)
  )
}

### Graphical exploratory data analysis

#' Print one qqplot for each variable and group
#' 
#' @description `qqplot_by_factor` takes a data frame and prints a qq-plot for each group and varible
#' @details Take a data frame and prints a qq-plot for each group and varible. 
#' @param input_df a data frame with at least one factor column. 
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
  numeric_df <- select_if(input_df, is.numeric)
  
  func <- function(){
    par(mfrow=c(3,3))
    for (j in 1:ncol(numeric_df)){
      col_name <- colnames(numeric_df)[j]
      new_df <- na.omit(numeric_df)[,j]
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
      col_name <- colnames(numeric_df)[i]
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
#' @description `shapiro_by_factor` takes a data frame, applies the shapiro-wilk test 
#' to every group and variable and returns a table of p-values. 
#' @param input_df a data frame with at least one factor column
#' @param factor a string with the column name to group by
#' @example 
#' shapiro_by_factor(iris, iris$Species)
shapiro_by_factor <- function(input_df, factor){

  shapiro_pvalue <- aggregate(dplyr::select_if(input_df, is.numeric), 
                              by = list(factor),
                              FUN = function(x) {y <- shapiro.test(x); c(y$p.value)})
  
  
  print("p < 0.05 ... no normal distribution; p > 0.05 ... normal distribution")
  shapiro_pvalue
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

#' Matrix heatmap 
#' 
#' @description `matrix_heatmap` displays values of a matrix in a heatmap
#' @param input_matrix matrix. For exaple correlation or matrix of ratios or differences.  
#' @param interactive logical. Print heatmap to device (FALSE, default) 
#' or open interactive heatmap in browser (TRUE). 
#' @param out_path string. Path to save heatmap to png. 
#' If out_path is empty the heatmap is printed to the device. 
#' @example 
#' iris_cor <- cor(iris[,-5], method = "spearman")
#' matrix_heatmap(iris_cor)
#' \dontrun{
#' dir.create(paste(getwd(), "/examples", sep = ""), showWarnings = FALSE)
#' dir <- paste(getwd(), "/examples/iris", sep = "")
#' matrix_heatmap(iris_cor, interactive = TRUE)
#' matrix_heatmap(iris_cor, out_path = dir)
#' }
matrix_heatmap <- function(input_df, 
                           title = "",
                           interactive = FALSE, 
                           out_path = "none"){
  melted_matrix <- reshape::melt(input_df)
  names(melted_matrix) <- c("x", "y", "value")
  head(melted_matrix)
  
  matrix_heatmap <- ggplot(data = melted_matrix, aes(x=x, y=y, fill=value)) +
    geom_tile() +
    ggtitle(title) +
    scale_fill_viridis_c(direction = -1) +
    my_theme +
    theme(axis.title.x = element_blank(), 
          axis.title.y = element_blank(),
          axis.text.y = element_text(size = 6),
          axis.text.x = element_text(angle = 90, size = 6, hjust = 1))
  
  if(out_path != "none"){
    print(paste("Saving plot to ", out_path, "_matrix_heatmap.png", sep = ""))
    ggsave(paste(out_path, "_matrix_heatmap.png", sep = ""), 
           plot = matrix_heatmap)
  }
  
  if (interactive == TRUE) {
    interactive_heatmap <- ggplotly(matrix_heatmap) # interactive heatmap
    interactive_heatmap
  } else {
    matrix_heatmap # static heatmap
  }
}

#' Calculate ratio matrix
#' 
#' @description `calculate_ratio_matrix` takes a vector and calculates the ratio between all elements of a vector. 
#' @param input_vector numeric vector.
#' @param names_vector vector. Length = length(input_vector). Row and column names for the returnded matrix. Default are numbers 1:length(input_vector). 
#' @example 
#' myvector <- c(1712.9583, 1446.4583, 1968.4167, 2124.1250, 2315.7083, 1135.1667, 2227.5000)
#' mynames <- letters[1:length(myvector)]
#' calculate_ratio_matrix(myvector, mynames)
calculate_ratio_matrix <- function(input_vector, 
                                   names_vector = 1:length(input_vector)){
  
  ratio_matrix <- matrix(nrow = length(input_vector), 
                         ncol = length(input_vector))
  
  for(i in 1:length(input_vector)){
    for(j in 1:length(input_vector)){
      buffer <- input_vector[i] / input_vector[j]
      ratio_matrix[i, j] <- buffer
    }
  }
  
  rownames(ratio_matrix) <- names_vector
  colnames(ratio_matrix) <- names_vector
  
  ratio_matrix
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
#' parallel_plot(iris, iris$Species)
#' parallel_plot(iris, 5, scale = "center", titles = c("Centered parallel plot", "Parameters", "Univariate scale to standardize vertical height"))
#' \dontrun
#' dir.create(paste(getwd(), "/examples", sep = ""), showWarnings = FALSE)
#' dir <- paste(getwd(), "/examples/iris", sep = "")
#' parallel_plot(iris, 5, out_path = dir)
parallel_plot <- function(input_df,  group_vector, out_path = "none", 
                          title = "Parallel Plot", 
                          x_title = "", 
                          y_title ="", 
                          legend_title = "", 
                          col_names = colnames(dplyr::select_if(input_df, is.numeric)), 
                          scale = "globalminmax"){

  new_df <- dplyr::select_if(input_df, is.numeric)
  new_df <- cbind(new_df, group_vector)
  
  par_plot <- ggparcoord(new_df,
                         columns = 1:(ncol(new_df)-1),
                         groupColumn = ncol(new_df),
                         showPoints = FALSE,
                         scale = scale, # "globalminmax" is default
                         alphaLines = 0.5
  ) +
    ggtitle(title) + 
    xlab(x_title) + 
    ylab(y_title) +
    labs(color = legend_title) +
    scale_color_viridis(discrete=TRUE, end = 1) +
    geom_point(shape = 20, size = 0.5) +
    theme(plot.title = element_text(size=12, hjust = 0.5, family="AvantGarde"),
          axis.text.x = element_text(angle = 45, size = 7, hjust = 1, colour = "grey40", family="AvantGarde"), 
          axis.text.y = element_text(size = 7, colour = "grey40", family="AvantGarde"),
          axis.title.x = element_text(size = 10, hjust = 0.5, colour = "grey40", family="AvantGarde"),
          axis.title.y = element_text(size = 10, hjust = 0.5, colour = "grey40", family="AvantGarde"),
          legend.text = element_text(size = 7, colour = "grey40", family="AvantGarde"),
          legend.title = element_text(size = 8, colour = "grey40", family="AvantGarde"), 
          legend.position = "right") 


  if(out_path != "none"){
    print(paste("Saving plot to ", out_path, "_parcoord.png", sep = ""))
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
#' spider_chart(minimized_iris, title = "Spider chart of iris", legend_lab = c(1, 2, 3))
#' \dontrun
#' minimized_iris <- aggregate(dplyr::select_if(iris, is.numeric), 
#'                             by = list(iris$Species), 
#'                             FUN = mean)
#' rownames(minimized_iris) <- minimized_iris$Group.1
#' dir.create(paste(getwd(), "/examples", sep = ""), showWarnings = FALSE)
#' dir <- paste(getwd(), "/examples/iris", sep = "")
#' spider_chart(minimized_iris, out_path = dir)
spider_chart <- function(minimized_df, 
                         title = "Spider chart", 
                         legend_lab = rownames(minimized_df),
                         out_path = "none"){ # todo get labels ot of the plot
  
  out_name <- paste(out_path, "_spiderChart", ".png", sep = "")
  
  func <- function(){
    par(mfrow = c(1, 1))
    spider_data <- dplyr::select_if(minimized_df, is.numeric) # remove column with rownames
    spider_labels <- colnames(spider_data) # set max. label length to 10 characters
    
    spider_min <- floor(min(spider_data))
    spider_max <- ceiling(max(spider_data))
    spider_interval <- (spider_max-spider_min)/4
    spider_data <- as.data.frame(select_if(spider_data, is.numeric))
    
    # add max and min to the dataframe to plot the grid
    spider_data <- rbind(spider_min, spider_max, spider_data)
    
    ## set colors
    colors_border = as.vector(viridis(n = nrow(minimized_df), option = "viridis"))
    colors_in = alpha(colors_border, alpha = 0.1)
    
    ## radar chart
    par(mfrow = c(1, 1))
    radarchart(spider_data,
               axistype=1,
               caxislabels = c(spider_min, 
                               spider_min+spider_interval, 
                               spider_min+2*spider_interval, 
                               spider_min+3*spider_interval, 
                               spider_min+4*spider_interval),
               #custom polygon
               pcol=colors_border,
               pfcol=colors_in,
               pty = 20,
               plwd=2,
               plty=1,
               # custom grid
               cglcol="grey",
               cglty=1,
               axislabcol="grey",
               cglwd=0.8,
               # custom labels
               vlcex=0.6, 
               vlabels = spider_labels,
               centerzero = FALSE)
    
    title(main = title, cex.main = 0.9, font.main = 1)
    
    ## Add a legend
    legend(x=-1.4, 
           y=1.1, 
           legend = legend_lab, 
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

#' Simple barchart
#' 
#' `simple_barchart` returns a ggplot barchart with viridis coloring
#' @param data_frame data frame with at least one numeric and one factor column
#' @param x vector. Factor column
#' @param y vector. Numeric column
#' @param fill vector. Factor column to fill bars by. By default all bars are blue. 
#' @param title string. Main title of the plot. Default = ""
#' @param xlab string. Title of x-axis. Default = ""
#' @param ylab string. Title of y-axis. Default = ""
#' @example 
#' iris_table <- reshape::melt(table(iris$Species))
#' simple_barchart(iris_table, x = iris_table$Var.1, y = iris_table$value, fill = "#35608DFF")
simple_barchart <- function(data_frame, x, y, fill = "#35608DFF", 
                            title = "test", xlab = "", ylab = ""){
  simple_barchart <- ggplot(data_frame, aes(x=x, y=y, fill=fill)) +
    geom_bar(stat="identity") + 
    labs(title = title, 
         x = xlab, 
         y = ylab) +
    scale_fill_viridis_d() +
    theme(legend.position = "none")
  
  
  if(is.character(fill)){
    simple_barchart <- simple_barchart + geom_bar(stat="identity", fill = fill)
  }
  simple_barchart
}

#' Plot ratio barcharts
#' @description `plot_ratio_barcharts` plots the ratios of a list of values from a data frame by group. 
#' @param data_frame data frame. From aggregate() function. 
#' @param subset_vector vector of column names from the data_frame. 
#' @param grop_column string. Name of the column to group by. 
#' @example 
#' aggregated_iris <- aggregate(iris[-5], by = list(iris$Species), FUN = mean)
#' plot_ratio_barcharts(aggregated_iris, c("Sepal.Length", "Sepal.Width", "Petal.Length"), group_column = "Group.1")
plot_ratio_barcharts <- function(data_frame, subset_vector, group_column){
  mysubset <- subset(data_frame, select = subset_vector)
  myratios <- lapply(1:nrow(mysubset), 
                     function(i) calculate_ratio_matrix(as.numeric(mysubset[i, ]), 
                                                        names_vector = colnames(mysubset)))
  names(myratios) <- as.character(data_frame[[group_column]])
  
  myratios <- reshape::melt(myratios)
  myratios$ratio_name <- as.factor(paste(myratios$X1, "/", myratios$X2))
  myratios <- myratios[-(1:2)]
  
  my_ratio_list <- lapply(1:length(levels(myratios$ratio_name)), 
                          function(i) subset(myratios, myratios$ratio_name == levels(myratios$ratio_name)[i]))
  names(my_ratio_list) <- levels(myratios$ratio_name)
  
  ratio_barcharts <- lapply(1:length(my_ratio_list), 
                            function(i) simple_barchart(my_ratio_list[[i]], 
                                                        x = my_ratio_list[[i]]$L1, 
                                                        y = my_ratio_list[[i]]$value,
                                                        fill = as.factor(my_ratio_list[[i]]$L1), 
                                                        title = my_ratio_list[[i]]$ratio_name, 
                                                        xlab = "group", 
                                                        ylab = "ratio"))
  ratio_barcharts
}
