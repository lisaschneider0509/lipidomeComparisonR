#' Hierarchical clustering heatmap
#' 
#' @description `hclust_heatmap` takes a data frame and prints a heatmap from hierarchical clustering of the numeric variables
#' @details This function takes a data frame and transforms it into a matrix. 
#' The given matrix is then used to print a heatmap with hierarchical clustering. The 
#' @examples 
#' hclust_heatmap(USArrests)
#' hclust_heatmap(iris, "Species")
hclust_heatmap <- function(input_df, 
                           set_rownames = "none", 
                           out_path = "none"){
  
  input_matrix <- as.matrix(dplyr::select_if(input_df, is.numeric))
  out_name <- paste(out_path, "_hierarchical_heatmap", ".png", sep = "")
  
  
  if(set_rownames != "none"){
    rownames(input_matrix) <- input_df[[set_rownames]]
  }
  
  func <- function(){
    par(mfrow = c(1, 1), cex.main = 0.8, cex.lab = 0.8)
    gplots::heatmap.2(input_matrix, 
                      scale = "row", # = c("none","row", "column")
                      
                      ## general apperance
                      trace = "none",
                      density.info = "none",
                      col = viridis::viridis_pal(), 
                      margins = c(5, 8),
                      
                      ## dendrogram and labels
                      cexRow = 0.6, 
                      cexCol = 0.7,
                      offsetRow = 0, 
                      offsetCol = 0, 
                      
                      ## other labels
                      main = mytitle)
  }
  
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
