## code to prepare `DATASET` dataset goes here

tfData <- read.csv('data-raw/tracefinder_example.csv')
usethis::use_data(tfData, overwrite = TRUE)
