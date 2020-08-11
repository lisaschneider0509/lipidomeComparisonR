#### R script with functions for lipidome comparison ####

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
#' @export
#' @example
#' calc_by_replicate(iris, iris$Species, mean)
calc_by_replicate <- function(input_df,
                              factor,
                              funct,
                              na_action = na.omit
){
  as.data.frame(
    stats::aggregate(dplyr::select_if(input_df, is.numeric),
              by=list(factor),
              FUN=funct,
              na.action = na_action)
  )
}

#' Generate table of categorical variable
#'
#' @description `generate_categorical_table` takes a categorical column of a data frame and generates a table where the categorical strings are assigned to a numeric code. This way the information of additional catecorical variabes is not lost when performing an aggregate function.
#' @param categorical_col column of data frame containing categorical values (i.e. strings). Format = df$colname
#' @export
#' @example
#' generate_categorical_table(iris$Species)
generate_categorical_table <- function(categorical_col){
  categorical_table <- cbind.data.frame(as.vector(categorical_col), as.numeric(categorical_col))
  categorical_table <- as.data.frame(table(categorical_table))
  categorical_table <- categorical_table[categorical_table$Freq != 0,][,-3]
  colnames(categorical_table) <- c("V1", "V2")
  categorical_table

}
