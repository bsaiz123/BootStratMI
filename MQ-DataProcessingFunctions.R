library("tibble")
library("dplyr")
library("tidyr")
library("purrr")
library("qvalue")
library("ggplot2")
library("stringr")
library("minet")
library("doSNOW")
library("foreach")
library("parallel")

# Data Prep
FilterContaminants <- function(df) {
  check_cols <- c("Reverse", "Only.identified.by.site", "Potential.contaminant")
  valid_cols <- check_cols[check_cols %in% colnames(df)]
  
  # If none of the specified columns exist, return the original dataframe
  if(length(valid_cols) == 0) {
    warning("None of the specified columns found in the dataframe")
    return(df)
  }

  keep_rows <- rep(TRUE, nrow(df))
  for(col in valid_cols) {
    keep_rows <- keep_rows & (df[[col]] != "+")
  }
  
  filtered_df <- df[keep_rows, ]
  
  removed_count <- nrow(df) - nrow(filtered_df)
  message(paste("Removed", removed_count, "rows containing '+' in the specified columns"))
  
  return(filtered_df)
}

LogTransformData <- function(df, dataCols, majorityProteinCol = "Majority.protein.IDs"){
  dfLog <- log2(df[,dataCols])
  dfLog[dfLog == -Inf] <- NA
  output <- cbind(df[, majorityProteinCol, drop = FALSE], dfLog)
  return(output)
}

FlagNA <- function(df, threshold = 0.7){
  # Return Boolean vector indicating which rows are at least threshold% complete
  row_completeness <- rowSums(!is.na(df)) / ncol(df)
  return(row_completeness >= threshold)
}

DefineTreatment <- function(df, treatments = c()){
  # Groups treatments by column names, uses identifiable strings
  treatList = list()
  for(i in 1:length(treatments)){
    treatList[[treatments[i]]] <- df[,grep(treatments[i], colnames(df), value = TRUE)]
  }
  return(treatList)
}

FilterFlagged <- function(df,
                          treatments = c(),
                          threshold = 0.7){
  # First split the data by treatments
  dataList <- DefineTreatment(df, treatments = treatments)
  
  # Create a vector to track which rows should be kept
  rowsToKeep <- rep(FALSE, nrow(df))
  
  # Check each treatment group
  for(treatment in names(dataList)){
    # Apply FlagNA to each treatment dataset
    validRows <- FlagNA(dataList[[treatment]], threshold = threshold)
    # Mark rows that are valid in any treatment group
    rowsToKeep <- rowsToKeep | validRows
  }
  
  # Keep only rows where at least one treatment group meets the threshold
  slicedData <- df[rowsToKeep, ]
  slicedData <- slicedData[, colSums(!is.na(slicedData)) > 0] # Remove empty columns
  
  return(slicedData)
}

ImputeMissingValues <- function(df, width = 0.3, downshift = 1.8, seed = 100){
  data <- as.data.frame(df)

  set.seed(seed)
  
  for (col in 1:ncol(data)) {
    # Skip non-numeric columns
    if (!is.numeric(data[[col]])) {
      message(paste("Skipping non-numeric column:", colnames(data)[col]))
      next
    }
    
    # Ensure the column is a simple numeric vector
    data[[col]] <- as.numeric(unlist(data[[col]]))
    
    # Get observed values in this column
    observed_values <- data[[col]][!is.na(data[[col]])]
    
    # Skip if no observed values
    if (length(observed_values) == 0) {
      message(paste("Skipping column with no observed values:", colnames(data)[col]))
      next
    }
    
    # Calculate mean and standard deviation of observed values
    mu <- mean(observed_values, na.rm = TRUE)
    sigma <- sd(observed_values, na.rm = TRUE)
    
    # Handle case where standard deviation is 0 or NA
    if (is.na(sigma) || sigma == 0) {
      message(paste("Column", colnames(data)[col], "has zero/NA standard deviation. Using small value."))
      sigma <- 0.01
    }
    
    # Calculate parameters for imputation distribution
    impute_mu <- mu - (downshift * sigma)
    impute_sigma <- width * sigma
    
    # Find which rows have missing values in this column
    missing_idx <- which(is.na(data[[col]]))
    
    # Generate random values from normal distribution
    if (length(missing_idx) > 0) {
      imputed_values <- rnorm(
        n = length(missing_idx),
        mean = impute_mu,
        sd = impute_sigma
      )
      
      # Replace missing values with imputed values
      data[missing_idx, col] <- imputed_values
    }
  }
  
  # Convert back to tibble
  result <- tibble::as_tibble(data)
  return(result)
}

RawDataPreparation <- function(df,
                               treatments,
                               dataColPattern = "^LFQ.intensity",
                               majorityProteinCol = "Majority.protein.IDs",
                               NA.Threshold = 0.7,
                               imputeWidth = 0.3,
                               imputeDownshift = 1.8,
                               seed = 100){
  
  df.filtered <- FilterContaminants(df)
  
  dataCols <- grep(dataColPattern, colnames(df), value = FALSE)
  
  df.LogTransformed <- LogTransformData(df.filtered, 
                                        dataCols = dataCols, 
                                        majorityProteinCol = majorityProteinCol)
  
  df.NA_Filtered <- FilterFlagged(df.LogTransformed, 
                                  treatments = treatments, 
                                  threshold = NA.Threshold)
  
  df.Imputed <- ImputeMissingValues(df.NA_Filtered,
                                    width = imputeWidth,
                                    downshift = imputeDownshift,
                                    seed = seed)
  return(df.Imputed)
}
