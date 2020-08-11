#' Print one qqplot for each variable and group
#'
#' @description `qqplot_by_factor` takes a data frame and prints a qq-plot for each group and varible
#' @details Take a data frame and prints a qq-plot for each group and varible.
#' @param input_df a data frame with at least one factor column.
#' @param factor a string with the column name to group by
#' @param out_path optional string. If out path is given, the generated plots are saved to a pdf document with the given name.
#' @export
#' @examples
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
    graphics::par(mfrow=c(3,3))
    for (j in 1:ncol(numeric_df)){
      col_name <- colnames(numeric_df)[j]
      new_df <- stats::na.omit(numeric_df)[,j]
      stats::qqnorm(numeric_df[,j][input_df[[factor]] == levels[i]],
             main = paste(col_name, levels[i], sep = " "),
             cex.main = 0.8)
      stats::qqline(numeric_df[,j][input_df[[factor]] == levels[i]])
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
      grDevices::pdf(paste(out_path, "_qqplot_", levels[i], ".pdf", sep = ""))
      func()
      grDevices::dev.off()
    }
  }
}

#' Print one histogram for each variable group.
#'
#' @description `histogram_by_factor` prints a histogram with density line for each group and variable of a data frame
#' @details Generate a new grid/pdf dochument for each group of a data frame. Generate one plot for each variable of this subset data frame.
#' @param input_df a data frame with at least one factor column
#' @param factor a string with the column name to group by
#' @param out_path optional string.
#' @export
#' @examples
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
    graphics::par(mfrow=c(3,3))
    for (j in 1:ncol(numeric_df)){
      col_name <- colnames(numeric_df)[j]
      graphics::hist(numeric_df[,j][input_df[[factor]] == levels[i]],
           main = paste(col_name, levels[i], sep = " "),
           cex.main = 0.8,
           xlab = NULL)
      graphics::lines(stats::density(numeric_df[,j][input_df[[factor]] == levels[i]]))
      graphics::lines(stats::density(numeric_df[,j][input_df[[factor]] == levels[i]],
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
      grDevices::pdf(paste(out_path, "_hist_", levels[i], ".pdf", sep = ""))
      func()
      grDevices::dev.off()
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
#' @export
#' @examples
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
      graphics::boxplot(numeric_df[,i] ~ input_df[[factor]],
              main = col_name,
              xlab = NULL,
              ylab = NULL)
    }
  }

  if(out_path == "none"){
    graphics::par(mfrow=c(3,3),
        cex.main = 1,
        cex.axis = 0.7)
    func()
  }
  else{
    print(paste("Saving to ", out_path, "_boxplot.pdf", sep = ""))
    grDevices::pdf(paste(out_path, "_boxplot.pdf", sep = ""))
    graphics::par(mfrow=c(3,3),
                  cex.main = 1,
                  cex.axis = 0.7)
    func()
    grDevices::dev.off()
  }
}


#' Assess normality for each group
#'
#' @description `shapiro_by_factor` takes a data frame, applies the shapiro-wilk test to every group and variable and returns a table of p-values.
#' @param input_df a data frame with at least one factor column
#' @param factor a string with the column name to group by
#' @export
#' @example
#' shapiro_by_factor(iris, iris$Species)
shapiro_by_factor <- function(input_df, factor){

  shapiro_pvalue <- stats::aggregate(dplyr::select_if(input_df, is.numeric),
                              by = list(factor),
                              FUN = function(x) {y <- stats::shapiro.test(x); c(y$p.value)})


  print("p < 0.05 ... no normal distribution; p > 0.05 ... normal distribution")
  shapiro_pvalue
}


#' Correlation plots
#'
#' @description `correlation_plot` calculates correlations of variables and displays them in a dotplot
#' @details A grid of dotplots displaying the correlations between the variables of a data frame is generated. Use with less than 10 variables ideally.
#' @param input_df data frame.
#' @param method string. Method for calculating the correlation. Options: "pearson", "kendall", "spearman" (default).
#' @param out_path optional string. Path to save correlation plot to png. If out_path is empty the correlation plot is printed to the device.
#' @export
#' @examples
#' correlation_plot(iris)
#' \dontrun{
#' dir.create(paste(getwd(), "/examples", sep = ""), showWarnings = FALSE)
#' dir <- paste(getwd(), "/examples/iris", sep = "")
#' correlation_plot(iris, out_path = dir)
#' }
correlation_plot <- function(input_df, method = "spearman", out_path = "none"){
  out_name <- paste(out_path, "_correlations", ".png", sep = "")

  func <- function(){
    graphics::par(cex.main = 1, cex.axis = 0.8, family = "AvantGarde", font = 1)
    colour = viridis::viridis(n = 2)
    psych::pairs.panels(dplyr::select_if(input_df, is.numeric),
                        method = method,
                        hist.col = colour[2], rug = FALSE,
                        density = TRUE,  lm = TRUE, ci = FALSE, col = colour[1],
                        ellipses = FALSE,
                        pch = 20, cex = 0.7,
                        cex.cor = 0.5, main = "Correlation plot")
  }

  if(out_path != "none"){
    print(paste("Saving to ", out_name, sep = ""))
    grDevices::png(out_name)
    func()
    grDevices::dev.off()
    func()
  }
  else{
    func()
  }
}

#' Matrix heatmap
#'
#' @description `matrix_heatmap` displays values of a matrix in a heatmap
#' @param input_df data frame or matrix. For exaple correlation or matrix of ratios or differences.
#' @param interactive logical. Print heatmap to device (FALSE, default) or open interactive heatmap in browser (TRUE).
#' @param title string. Main title of the plot. Default: no title
#' @param out_path string. Path to save heatmap to png. If out_path is empty the heatmap is printed to the device.
#' @export
#' @examples
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
  utils::head(melted_matrix)

  matrix_heatmap <- ggplot2::ggplot(data = melted_matrix,
                                    ggplot2::aes(x=melted_matrix$x,
                                                 y=melted_matrix$y,
                                                 fill=melted_matrix$value)) +
    ggplot2::geom_tile() +
    ggplot2::ggtitle(title) +
    ggplot2::scale_fill_viridis_c(direction = -1) +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(size=12, hjust = 0.5, family="AvantGarde"),
                   axis.text.x = ggplot2::element_text(angle = 90, size = 6, hjust = 1, colour = "grey40", family="AvantGarde"),
                   axis.text.y = ggplot2::element_text(size = 6, colour = "grey40", family="AvantGarde"),
                   axis.title.x = ggplot2::element_blank(),
                   axis.title.y = ggplot2::element_blank()
    )

  if(out_path != "none"){
    print(paste("Saving plot to ", out_path, "_matrix_heatmap.png", sep = ""))
    ggplot2::ggsave(paste(out_path, "_matrix_heatmap.png", sep = ""),
           plot = matrix_heatmap)
  }

  if (interactive == TRUE) {
    interactive_heatmap <- plotly::ggplotly(matrix_heatmap) # interactive heatmap
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
#' @export
#' @examples
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
#' @details This function takes a data frame with at least one factor variable and displays it in a pralell coordinates plot, where the different groups are color coded. This plot works best for <= 10 parameters.
#' @param input_df data frame.
#' @param group_vector vector. Vector (part of the data frame), to sort by.
#' @param out_path string. Path to save parallel plot to png. If out_path is empty the parallel plot is printed to the device.
#' @param title string. Main Title. Default = "Parallel plot"
#' @param x_title string. Title of x-axis. Default: none
#' @param y_title string. Title of y-axis. Default: none
#' @param legend_title string. Title of legend. Default: none
#' @param col_names vector. Column names of the data frame and lables on the x-axis. Default: colnames(input_df)
#' @param scale string. Method used to scale the variable. Default = "globalminmax". Options: "std", "robust", "uniminmax", "globalminmax", "center", "centerObs". For more information on the options see help(ggparcoord).
#' @export
#' @examples
#' parallel_plot(iris, iris$Species)
#' \dontrun{
#' dir.create(paste(getwd(), "/examples", sep = ""), showWarnings = FALSE)
#' dir <- paste(getwd(), "/examples/iris", sep = "")
#' parallel_plot(iris, 5, out_path = dir)}
parallel_plot <- function(input_df,
                          group_vector,
                          out_path = "none",
                          title = "Parallel Plot",
                          x_title = "",
                          y_title ="",
                          legend_title = "",
                          col_names = colnames(dplyr::select_if(input_df, is.numeric)),
                          scale = "globalminmax"){

  new_df <- dplyr::select_if(input_df, is.numeric)
  new_df <- cbind(new_df, group_vector)

  par_plot <- GGally::ggparcoord(new_df,
                         columns = 1:(ncol(new_df)-1),
                         groupColumn = ncol(new_df),
                         showPoints = FALSE,
                         scale = scale, # "globalminmax" is default
                         alphaLines = 0.5
  ) +
    ggplot2::ggtitle(title) +
    ggplot2::xlab(x_title) +
    ggplot2::ylab(y_title) +
    ggplot2::labs(color = legend_title) +
    viridis::scale_color_viridis(discrete=TRUE, end = 1) +
    ggplot2::geom_point(shape = 20, size = 0.5) +
    ggplot2::theme(plot.title = ggplot2::element_text(size=12, hjust = 0.5, family="AvantGarde"),
          axis.text.x = ggplot2::element_text(angle = 45, size = 7, hjust = 1, colour = "grey40", family="AvantGarde"),
          axis.text.y = ggplot2::element_text(size = 7, colour = "grey40", family="AvantGarde"),
          axis.title.x = ggplot2::element_text(size = 10, hjust = 0.5, colour = "grey40", family="AvantGarde"),
          axis.title.y = ggplot2::element_text(size = 10, hjust = 0.5, colour = "grey40", family="AvantGarde"),
          legend.text = ggplot2::element_text(size = 7, colour = "grey40", family="AvantGarde"),
          legend.title = ggplot2::element_text(size = 8, colour = "grey40", family="AvantGarde"),
          legend.position = "right") +
    ggplot2::theme_minimal()


  if(out_path != "none"){
    print(paste("Saving plot to ", out_path, "_parcoord.png", sep = ""))
    ggplot2::ggsave(paste(out_path, "_parcoord.png", sep = ""),
           plot = par_plot)
  }
  else{
    par_plot
  }
}


#' Spider Chart
#'
#' @description `spider chart` takes a minimized data frame (one value per group) and prints a spider chart
#' @details {This function takes a data frame of with one value per group (i.e. calculate mean groupwise).
#' This minimized data frame is used to draw a spider chart (also radar chart odr network plot).
#' The ideal numer of parameters for a spider chart is <= 10. Also the shape of the graph depends on
#' the order of parameters. If the data has large differences in size, normalizing or scaling the data is necessary.
#' Doesn't work with non normal data.}
#' @param minimized_df data frame. A data frame that contains only one value per group and variable. Most often this will be a data frame of means calculated from another data frame.
#' @param title string. Main title of the chart. Default = "Spider chart"
#' @param legend_lab vector. Lables. Default: rownames(minimized_df)
#' @param out_path string. Path to save spider chart to png. If out_path is empty, the spider chart is printed to the device.
#' @export
#' @examples
#' minimized_iris <- aggregate(dplyr::select_if(iris, is.numeric),
#'                             by = list(iris$Species),
#'                             FUN = mean)
#' rownames(minimized_iris) <- minimized_iris$Group.1
#' spider_chart(minimized_iris, title = "Spider chart of iris", legend_lab = c(1, 2, 3))
spider_chart <- function(minimized_df,
                         title = "Spider chart",
                         legend_lab = rownames(minimized_df),
                         out_path = "none"){

  out_name <- paste(out_path, "_spiderChart", ".png", sep = "")

  func <- function(){
    graphics::par(mfrow = c(1, 1))
    spider_data <- dplyr::select_if(minimized_df, is.numeric) # remove column with rownames
    spider_labels <- colnames(spider_data) # set max. label length to 10 characters

    spider_min <- floor(min(spider_data))
    spider_max <- ceiling(max(spider_data))
    spider_interval <- (spider_max-spider_min)/4
    spider_data <- as.data.frame(dplyr::select_if(spider_data, is.numeric))


    spider_data <- rbind(spider_min, spider_max, spider_data) # add max and min to the dataframe to plot the grid


    colors_border = as.vector(viridis::viridis(n = nrow(minimized_df), option = "viridis")) ## set colors
    colors_in = scales::alpha(colors_border, alpha = 0.1)


    graphics::par(mfrow = c(1, 1)) # radar chart
    fmsb::radarchart(spider_data,
               axistype=1,
               caxislabels = c(spider_min,
                               spider_min+spider_interval,
                               spider_min+2*spider_interval,
                               spider_min+3*spider_interval,
                               spider_min+4*spider_interval),

               pcol=colors_border, #custom polygon
               pfcol=colors_in,
               pty = 20,
               plwd=2,
               plty=1,

               cglcol="grey", # custom grid
               cglty=1,
               axislabcol="grey",
               cglwd=0.8,
               # custom labels
               vlcex=0.6,
               vlabels = spider_labels,
               centerzero = FALSE)

    graphics::title(main = title, cex.main = 0.9, font.main = 1)


    graphics::legend(x=-1.4, # Add a legend
           y=1.1,
           legend = legend_lab,
           bty = "n",
           pch=20,
           col=colors_border,
           text.col = "black",
           cex=0.7, pt.cex=1.3)
  }


  if(out_path != "none"){ # print to device or save
    print(paste("Saving to ", out_name))
    grDevices::png(filename = out_name)
    func()
    grDevices::dev.off()
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
  simple_barchart <- ggplot2::ggplot(data_frame, ggplot2::aes(x=x, y=y, fill=fill)) +
    ggplot2::geom_bar(stat="identity") +
    ggplot2::labs(title = title,
         x = xlab,
         y = ylab) +
    viridis::scale_fill_viridis_d() +
    ggplot2::theme(plot.title = ggplot2::element_text(size=12, hjust = 0.5, family="AvantGarde"),
                   axis.text.x = ggplot2::element_text(size = 8, colour = "grey40", family="AvantGarde"),
                   axis.text.y = ggplot2::element_text(size = 8, colour = "grey40", family="AvantGarde"),
                   axis.title.x = ggplot2::element_text(size = 10, hjust = 0.5, colour = "grey40", family="AvantGarde"),
                   axis.title.y = ggplot2::element_text(size = 10, hjust = 0.5, colour = "grey40", family="AvantGarde"),
                   legend.position = "none") +
    ggplot2::theme_minimal()


  if(is.character(fill)){
    simple_barchart <- simple_barchart + ggplot2::geom_bar(stat="identity", fill = fill)
  }
  simple_barchart
}

#' Plot ratio barcharts
#' @description `plot_ratio_barcharts` plots the ratios of a list of values from a data frame by group.
#' @param data_frame data frame. From aggregate() function.
#' @param subset_vector vector of column names from the data_frame.
#' @param group_column string. Name of the column to group by.
#' @example
#' aggregated_iris <- aggregate(iris[-5], by = list(iris$Species), FUN = mean)
#' plot_ratio_barcharts(aggregated_iris, c("Sepal.Length", "Sepal.Width", "Petal.Length"), group_column = "Group.1")
plot_ratio_barcharts <- function(data_frame,
                                 subset_vector,
                                 group_column){
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
