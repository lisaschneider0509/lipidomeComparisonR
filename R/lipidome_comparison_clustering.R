#' Hierarchical clustering performance table
#'
#' @description `hclust_performance_table` takes a data frame and calculates which distance function and hclust function combination has the best performance.
#' @details This function calculated the performance for all combinations of the following distance and hclust functions distance functions = c("euclidean", "manhattan"), hclust functions = c("average", "single", "complete")). The results are printed in a table. The combination with the highest performance value will produce the best clustering.
#' @param input_df data frame.
#' @param out_path string. Path to save the text output to. If not specified, the output is only printed to the device or saved to a variable. If specified, the output is additionally saved in a .txt-file.
#' @param dist_methods vector. Options: c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"). Default: c("euclidean", "maximum", "manhattan", "canberra", "minkowski")
#' @param hclust_methods vector. Options: c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid"). Default: c("average", "single", "complete")
#' @export
#' @example
#' hclust_performance_table(USArrests)
#' \dontrun{
#' dir.create(paste(getwd(), "/examples", sep = ""), showWarnings = FALSE)
#' dir <- paste(getwd(), "/examples/USArrests", sep = "")
#' hclust_performance_table(USArrests, out_path = dir)
#' }
hclust_performance_table <- function(input_df,
                                     dist_methods = c("euclidean", "maximum", "manhattan", "canberra", "minkowski"),
                                     hclust_methods = c("average", "single", "complete"),
                                     out_path = "none"){

  buffer <- dendextend::dend_expend(input_df,
                        dist_methods = dist_methods,
                        hclust_methods = hclust_methods)
  performance <- buffer$performance

  if(out_path != "none"){
    print(paste("Writing to", out_path, "_hclust_performance.txt", sep = ""))
    utils::write.table(performance, paste(out_path, "_hclust_performance.txt", sep = ""))
    grDevices::dev.off()
  }
  performance
}

#' Hierarchical clustering performance plot
#'
#' @description `hclust_performance_plot` takes a data frame and calculates which distance function and hclust function combination has the best performance.
#' @details This function calculated the performance for all combinations of the following distance and hclust functions distance functions = c("euclidean", "manhattan"), hclust functions = c("average", "single", "complete")). The results are printed in a table. The combination with the highest performance value will produce the best clustering.
#' @param input_df data frame.
#' @param out_path string. Path to save the text output to. If specified, the output is saved in a .txt-file. If not specified, the output is printed th the device.
#' @param dist_methods vector. Options: c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"). Default: c("euclidean", "maximum", "manhattan", "canberra", "minkowski")
#' @param hclust_methods vector. Options: c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid"). Default: c("average", "single", "complete")
#' @export
#' @example
#' hclust_performance_plot(USArrests)
#' \dontrun{
#' dir.create(paste(getwd(), "/examples", sep = ""), showWarnings = FALSE)
#' dir <- paste(getwd(), "/examples/USArrests", sep = "")
#' hclust_performance_plot(USArrests, out_path = dir)
#' }
hclust_performance_plot <- function(input_df,
                                    dist_methods = c("euclidean", "manhattan"),
                                    hclust_methods = c("average", "single", "complete"),
                                     out_path = "none"){

  numeric_df <- dplyr::select_if(input_df, is.numeric)

  buffer <- dendextend::dend_expend(numeric_df,
                        dist_methods = dist_methods,
                        hclust_methods = hclust_methods)

  performance <- buffer$performance
  performance <- performance[colSums(!is.na(performance)) > 0]


  myplot <- ggplot2::ggplot(data=performance, aes(x=hclust_methods, y=optim, fill = dist_methods)) +
    ggplot2::geom_bar(stat="identity", position=position_dodge()) +
    ggplot2::geom_hline(yintercept = max(performance$optim), colour = "grey40", size = 0.3, linetype = "dashed") +
    ggplot2::scale_fill_viridis(discrete = TRUE, name = "Distance function") +
    ggplot2::labs(title="Hierarchical clustering methods - performance",
         subtitle = "",
         x="Hierarchical clustering function", y = "Performance") +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = element_text(size=12, hjust = 0.5),
              axis.text.x = element_text(size = 8),
              axis.title.x = element_text(size = 10, hjust = 0.5),
              axis.title.y = element_text(size = 10, hjust = 0.5),
              legend.text = element_text(size = 8),
              legend.title = element_text(size = 10),
              legend.position = "bottom")

  myplot
}


#' Hierarchical clustering method comparison
#'
#' @description `hclust_methods` takes a data frame and prints three dendrograms, one each for average linkage, single linkage and complete linkage.
#' @details This function takes a data frame and transforms it into a matrix. The given matrix is then used to perform hiererchical clustering with the average, single and complete linkage methods. The three clustering methods are then displayed as dendrograms. This plot is ment as a help to find ozt which clustering method should be used.
#' @param input_df data frame.
#' @param title string. Plot title. Default = "Hierarchical clustering method comparison".
#' @param dist_method string. MEthod for distance function. = c("euclidean", "manhattan")
#' @param out_path string. Path to save plot to png. Default = to device.
#' @param labs vector. Labels of dendrogram. Default: rownames(input_df)
#' @export
#' @examples
#' hclust_methods(USArrests)
#' hclust_methods(USArrests, labs = 1:nrow(USArrests), dist_method = "manhattan")
#' \dontrun{
#' dir.create(paste(getwd(), "/examples", sep = ""), showWarnings = FALSE)
#' dir <- paste(getwd(), "/examples/USArrests", sep = "")
#' hclust_methods(USArrests, out_path = dir)
hclust_methods <- function(input_df,
                           labs = row.names(input_df),
                           title = "Hierarchical clustering method comparison",
                           dist_method = "euclidean",
                           out_path = "none"){

  out_name <- paste(out_path, "_hclust_methods.png", sep = "")
  average_hclust <- stats::hclust(stats::dist(dplyr::select_if(input_df, is.numeric), method = dist_method),
                           method = "average")
  single_hclust <- stats::hclust(stats::dist(dplyr::select_if(input_df, is.numeric), method = dist_method ),
                          method = "single")
  complete_hclust <- stats::hclust(stats::dist(dplyr::select_if(input_df, is.numeric), method = dist_method),
                            method = "complete")
  mycolors <- viridis::viridis(n = 4)

  func <- function(){
    graphics::par(mfrow=c(1,3), oma = c(0, 0, 3, 0), cex = 0.7, family = "AvantGarde")
    graphics::plot(average_hclust, main = "Average linkage",
         labels = labs,
         xlab = "", sub = "",
         col = mycolors[1])
    graphics::plot(single_hclust, main = "Single linkage",
         xlab = "", sub = "",
         labels = labs,
         col = mycolors[2])
    graphics::plot(complete_hclust, main = "Complete linkage",
         xlab = "", sub = "",
         labels = labs,
         col = mycolors[3])
    graphics::mtext(title,
          outer = TRUE,
          cex = 1,
          # side=3,
          line=1,
          family = "AvantGarde",
          font = 1
          )
    graphics::mtext(paste(dist_method, " distance", sep = ""),
          # side=1,
          outer = TRUE,
          line=0,
          family = "AvantGarde",
          font = 3, cex=0.8)

    # graphics::mtext(title, side=3, line=2, family = "AvantGarde", font = 1, cex = 1)
  }

  if(out_path != "none"){
    print(paste("Saving to ", out_name))
    grDevices::png(filename = out_name, width = 1920, height = 600, unit = "px")
    func()
    grDevices::dev.off()
  }
  else{
    func()
  }

}

#' Hierarchical clustering dendrogram
#'
#' @description `hclust_dendrogram` takes a data frame and prints a dendrogram.
#' @details This function takes a data frame and transforms it into a matrix. The given matrix is then used to perform hiererchical clustering with the average, single or complete linkage. The clustering is plotted in a dendrogram.
#' @param hclust_element object of class hclust.
#' @labs vector or list. Labels fr the branches of the dendrogram. Default are the labels from hclust_element$labels.
#' @param title string. Plot title. Default = "Hierarchical clustering".
#' @param out_path string. Path to save plot to png. Default = to device.
#' @param labs vector. Labels of dendrogram. Default: hclust_element$labels
#' @export
#' @examples
#' US_clust <- stats::hclust(stats::dist(USArrests, method = "euclidean"), method = "average")
#' hclust_dendrogram(US_clust)
#' hclust_dendrogram(US_clust, labs = 1:nrow(USArrests))
#' \dontrun{
#' dir.create(paste(getwd(), "/examples", sep = ""), showWarnings = FALSE)
#' dir <- paste(getwd(), "/examples/USArrests", sep = "")
#' hclust_dendrogram(US_clust, out_path = dir)}
hclust_dendrogram <- function(hclust_element,
                              labs = hclust_element$labels,
                              title = "Hierarchical clustering",
                              out_path = "none"){

  out_name <- paste(out_path, "_dendrogram_", hclust_element$method, "_linkage.png", sep = "")
  mycolors <- viridis::viridis(n = 4, begin = 0, end = 0)

  linkage <- hclust_element$method
  distance <- hclust_element$dist.method
  subtitle <- paste(distance, "distance, ", linkage, "linkage", sep = " ")

  func <- function(){
    graphics::par(mfrow = c(1, 1), cex = 1)
    graphics::plot(hclust_element, labels = labs,
         xlab = "",
         main = "", sub = "",
         col = "black")
    graphics::mtext(title, side=3, line=2, family = "AvantGarde", font = 1, cex = 1)
    graphics::mtext(subtitle, side=3, line=0.5, family = "AvantGarde", font = 3, cex=0.8)
  }

  if(out_path != "none"){
    print(paste("Saving to ", out_name))
    grDevices::png(filename = out_name, width = 1020, height = 600, unit = "px")
    func()
    grDevices::dev.off()
  }
  else{
    func()
  }
}

#' Hierarchical clustering heatmap
#'
#' @description `hclust_heatmap` takes a data frame and prints a heatmap from hierarchical clustering of the numeric variables
#' @details This function takes a data frame and transforms it into a matrix. The given matrix is then used to print a heatmap with hierarchical clustering.
#' @param input_df data frame.
#' @param row_names vector of length nrow(input_df). Default = rownames(input_df).
#' @param col_names vector of length ncol(select_if(input_df, is.numeric)). Default = colnames(select_if(input_df, is.numeric))
#' @param dist_method string. Method for the dist() function. dist_method = c("euclidean", "manhattan"). Default = "euclidean".
#' @param hclust_method string. Method for the hclust() function. hclust_method = c("complete", "single", "average"). Default = "complete".
#' @param title string. Plot title. Default = "Hierarchical clustering".
#' @param out_path string. Path to save plot to png.
#' @export
#' @examples
#' hclust_heatmap(USArrests)
#' hclust_heatmap(iris, row_names = iris$Species, dist_method = "manhattan", hclust_method = "average")
hclust_heatmap <- function(input_df,
                           row_names = rownames(input_df),
                           col_names = colnames(dplyr::select_if(input_df, is.numeric)),
                           dist_method = "euclidean",
                           hclust_method = "complete",
                           title = "Hierarchical clustering",
                           out_path = "none"){

  input_matrix <- as.matrix(dplyr::select_if(input_df, is.numeric))

  rownames(input_matrix) <- row_names
  colnames(input_matrix) <- col_names

  func <- function(){
    graphics::par(mfrow = c(1, 1),
        cex.main = 0.8, cex.lab = 0.8,
        family = "AvantGarde", font = 1)
    gplots::heatmap.2(input_matrix,
                      scale = "row", # = c("none","row", "column"),
                      # distfun = function(x) stats::dist(x, method="euclidean"),
                      # hclustfun = function(x) stats::hclust(x, method="ward.D2"),

                      ## general apperance
                      trace = "none",

                      # color key + density info
                      key = TRUE,
                      density.info = "none",
                      col = viridis::viridis_pal(),
                      margins = c(5, 8),

                      ## dendrogram and labels
                      cexRow = 0.6,
                      cexCol = 0.7,
                      offsetRow = 0,
                      offsetCol = 0,

                      ## other labels
                      main = title
                      )
  }

  if(out_path != "none"){
    print(paste("Saving to ", out_path))
    grDevices::png(filename = out_path, width = 16, height = 9, units = "in", res = 150)
    func()
    grDevices::dev.off()
  }
  else{
    func()
  }

}


#' Interactive hclust heatmap
#'
#' @description `hclust_heatmap_interactive` takes a data frame and lets you open an interactive heatmap in the viewer or save the heatmap to png.
#' @details This finction takes a data frame (optionally numeric or numeric with additional non-numeic varibles). Non-numeric variables are not part of the heatmap, but displayed in a separate box where the levels are color coded. The used method to create the heatmap is hierarchical clustering. Dendrograms are displayed for bowth rows and columns.
#' @param input_df data frame.
#' @param row_names data frame. Default are the row names of the input data frame.
#' @param col_names data frame. Default are the column names of the input data frame.
#' @param dist_method string. Method for the dist() function. dist_method = c("euclidean", "manhattan"). Default = "euclidean".
#' @param hclust_method string. Method for the hclust() function. hclust_method = c("complete", "single", "average"). Default = "complete".
#' @param title string. Plot title. Default = "Hierarchical clustering".
#' @param out_path string. Path to save plot to png.
#' @export
#' @examples
#' hclust_heatmap_interactive(iris, iris$Species)
#' hclust_heatmap_interactive(USArrests,
#'                            col_names = c("M", "A", "UP", "R"),
#'                            row_names = 1:50,
#'                            title = "Heatmap of USArrests")
#' \dontrun{
#' dir <- "examples/USArrests"
#' hclust_heatmap_interactive(USArrests, out_path = dir)}
hclust_heatmap_interactive <- function(input_df,
                                       row_names = rownames(input_df),
                                       col_names = colnames(dplyr::select_if(input_df, is.numeric)),
                                       dist_method = "euclidean",
                                       hclust_method = "complete",
                                       title = "Hierarchical clustering",
                                       # html_path = "",
                                       out_path = "none"
                                       ){

  out_file <- paste(out_path, "_interactive_hcheatmap.png", sep = "")

  hheatmap <- heatmaply::heatmaply(input_df,
                        ## clustering settings
                        dist_method = dist_method,
                        hclust_method = hclust_method,

                        ## display method
                        plot_method = "plotly", # = c("ggplot", "plotly"),

                        ## heatmap and dendrogram settings
                        dendrogram = "both", # = c("both", "row", "column", "none")
                        scale = "row", # = c("none","row", "column")
                        branches_lwd = 0.2,

                        ## colorbar & legend
                        hide_colorbar = FALSE,
                        colorbar_xpos = 1.02,
                        colorbar_ypos = 0.09,

                        ## Labels
                        main = title,
                        # sub = paste(dist_method, hclust_method, sep = "/"), #todo find out if and how subtitles work in heatmaply
                        xlab = "",
                        ylab = "",
                        labCol = col_names, # defaults to colnames; only takes the names of numeric colums as colnames
                        labRow = row_names, # defaults to rownames

                        ## general appearance
                        margins = c(60,100,40,20),
                        fontsize_row = 10, fontsize_col = 10,
                        heatmap_layers = theme(axis.line=element_blank()) #,
                        # file = html_path

  )

  if(out_path != "none"){
    print(paste("Saving to ", out_file, sep = ""))
    heatmaply::orca(hheatmap, file = out_file, width = 16 * 96, height = 8 * 96)
  }

  hheatmap
}
