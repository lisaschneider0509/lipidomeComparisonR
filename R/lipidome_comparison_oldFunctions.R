#' Takes a data frame with specific rownames and extracts meta data from it
#' 
#' @description 
#' `SID_to_metadata` takes a data frame with rownames in the form of 
#' <experiment_treatment_biol.rep.nr_tech.rep.nr> 
#' and extracts meta data
#' @details 
#' This function takes the information saved in sample IDs of the form 
#' <experiment_treatment_biol.rep.nr_tech.rep.nr> and saves each 
#' part of the sample ID as a column of character strings. 
#' If one part is missing, the column is filled with NAs. 
#' @param input_df a data frame with character strings and numerics
#' @examples 
#' df <- as.data.frame(rbind(c(1, 2, 3, 4), c(5, 6, 7, 8)))
#' row.names(df) <- c("Ex1_Tr1_1_1", "Ex1_Tr1_1_2")
#' colnames(df) <- c("a", "b", "c", "d")
#' ndf <- SID_to_metadata(df); ndf
SID_to_metadata <- function(input_df){
  n_sep <- mean(stringr::str_count(row.names(input_df), "_"))
  meta_info <- read.table(text = row.names(input_df), sep = "_")
  print(meta_info)
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


#' Generate table of categorical variable
#' 
#' @description `generate_categorical_table` takes a categorical column of a data frame and generates a table where the 
#' categorical strings are assigned to a numeric code. This way the information of additional catecorical variabes is not lost when performing 
#' an aggregate function. 
#' @param categorical_col column of data frame containing categorical values (i.e. strings). Format = df$colname
generate_categorical_table <- function(categorical_col){
  categorical_table <- cbind.data.frame(as.vector(categorical_col), as.numeric(categorical_col))
  categorical_table <- as.data.frame(table(categorical_table))
  categorical_table <- categorical_table[categorical_table$Freq != 0,][,-3]
  colnames(categorical_table) <- c("V1", "V2")
  categorical_table
  
}

#' Paste categroical variable
#' 
#' @description `paste_catecorical_variable` takes a dataframe produced by the aggregate function, the column number of a categorical variable, that was transformed
#' to numeric before applying aggregate and a table with the information of the numeric codes and the strings of the categorical variable and pastes the strings into 
#' the aggregated data frame. This function works only if the categorial variable has less different attributes than the variable used in the aggregate function. 
#' @param aggregated_df data frame produced by the aggregate() function.
#' @param categorical_col integer. Column of the aggregated_df that contains a categorical variable that was transformed to numeric before applying aggregate. 
#' @param categorical_table data frame. Contains information <which numeric value was assigned to which categorical attribute string. Produced by generate_categorical_table() function. 
#' 
paste_catecorical_variable <- function(aggregated_df, categrical_col, categorical_table){
  
  new_df <- data.frame()
  for(i in 1:nrow(categorical_table)){
    for(j in 1:nrow(aggregated_df)){
      if(categorical_table$V2[i] == aggregated_df[j, categrical_col]){
        buffer <- data.frame(aggregated_df$Group.1[j], categorical_table$V1[i])
        new_df <- rbind(new_df, buffer)
      }
    }
  }
  colnames(new_df) <- c("V1", "V2")
  new_df <- new_df[with(new_df, order(V1)),  ]
  aggregated_df <- aggregated_df[with(aggregated_df, order(Group)), ]
  aggregated_final <- aggregated_df
  aggregated_final[, categrical_col] <- new_df$V2
  aggregated_final
}

#' Normalize by median
med_normalize <- function(data_frame){
  norm_lipids <- select_if(data_frame, is.numeric)
  med_lipids <- apply(norm_lipids, 2, median, na.rm = TRUE)
  
  div <- function(x, y){
    z <- x / y
    z
  }
  
  new_df <- norm_lipids
  for(i in 1:ncol(norm_lipids)){
    buffer <- mapply(div, norm_lipids[i], med_lipids[i])
    new_df[,i] <- buffer
  }
  new_df
}
