# Make Metadata
ExtractMetadata <- function(df,
                            treatments,
                            timepoints,
                            ignoreText = "") {
  
  BuildPattern <- function(target) {
    if (length(ignoreText) == 2) {
      paste0(ignoreText[1], target, ignoreText[2])
    } else {
      paste0(ignoreText, target, "([^0-9]|$)")
    }
  }
  
  intensity_cols <- grep("^LFQ.intensity", colnames(df), value = TRUE)
  
  metadata <- data.frame(sample = intensity_cols, stringsAsFactors = FALSE)
  
  metadata$treatment <- sapply(intensity_cols, function(col) {
    match <- treatments[sapply(treatments, function(trt) grepl(BuildPattern(trt), col))]
    if (length(match) == 1) match else NA
  })
  
  metadata$timepoint <- sapply(intensity_cols, function(col) {
    match <- timepoints[sapply(timepoints, function(tp) grepl(BuildPattern(tp), col))]
    if (length(match) == 1) as.numeric(gsub("\\D", "", match)) else NA
  })
  
  unmatched_trt <- sum(is.na(metadata$treatment))
  unmatched_tp <- sum(is.na(metadata$timepoint))
  
  if (unmatched_trt > 0) warning(paste(unmatched_trt, "samples could not be matched to a treatment group"))
  if (unmatched_tp > 0) warning(paste(unmatched_tp, "samples could not be matched to a timepoint"))
  
  return(metadata)
}

# Calculate MI

TimepointCountData <- function(metadata, treatmentColumn = "treatment") {
  timepointSummary <- metadata %>%
    group_by(timepoint, .data[[treatmentColumn]]) %>%
    summarise(count = n(), .groups = 'drop') %>%
    pivot_wider(names_from = .data[[treatmentColumn]], 
                values_from = count, 
                names_prefix = "n_",
                values_fill = 0)
  return(timepointSummary)
}

StratifiedBootstrapIndices <- function(metadata, timepointCountData, treatmentColumn = "treatment") {

  treatments <- unique(metadata[[treatmentColumn]])
  sampled <- setNames(vector("list", length(treatments)), treatments)
  
  for(i in 1:nrow(timepointCountData)) {
    tp <- timepointCountData$timepoint[i]
    
    for(trt in treatments) {
      n <- timepointCountData[[paste0("n_", trt)]][i]
      tp_indices <- which(metadata$timepoint == tp & metadata[[treatmentColumn]] == trt)
      sampled[[trt]] <- c(sampled[[trt]], sample(tp_indices, size = n, replace = TRUE))
    }
  }
  return(sampled)
}

BootMI <- function(data,
                   metadata,
                   treatmentColumn, 
                   n_bootstrap,
                   treatments = c("DD","WW"),
                   estimator = "spearman",
                   parallel = FALSE,
                   workers = 0,
                   timeStrat = FALSE){
  
  data <- as.data.frame(data)
  rownames(data) <- data$Majority.protein.IDs
  data <- data[, sapply(data, is.numeric)]
  data <- t(data)
  
  group1_indices <- which(metadata[[treatmentColumn]] == treatments[1])
  group2_indices <- which(metadata[[treatmentColumn]] == treatments[2])
  
  timeStratifiedMetadata <- TimepointCountData(metadata,treatmentColumn)
  
  group1_data <- data[group1_indices, ]
  group2_data <- data[group2_indices, ]
  
  # Calculate observed matrices once
  OBS_mi_matrix1 <- build.mim(group1_data, estimator = estimator)
  OBS_mi_matrix2 <- build.mim(group2_data, estimator = estimator)
  
  # Count matrices for comparisons
  obs1_vs_boot2_higher <- matrix(0, nrow = ncol(data), ncol = ncol(data))
  obs2_vs_boot1_higher <- matrix(0, nrow = ncol(data), ncol = ncol(data))
  obs1_vs_boot2_lower <- matrix(0, nrow = ncol(data), ncol = ncol(data))
  obs2_vs_boot1_lower <- matrix(0, nrow = ncol(data), ncol = ncol(data))
  
  if(!parallel){
    StartTime <- Sys.time()
    for(i in 1:n_bootstrap){
      set.seed(i)
      if(!timeStrat){
        bootIndex1 <- sample(1:nrow(group1_data), size = nrow(group1_data), replace = TRUE)
        bootIndex2 <- sample(1:nrow(group2_data), size = nrow(group2_data), replace = TRUE)
        
        bootData1 <- group1_data[bootIndex1, ]
        bootData2 <- group2_data[bootIndex2, ]
      }else{
        bootstrapIndices <- StratifiedBootstrapIndices(metadata, timeStratifiedMetadata, treatmentColumn)
        bootData1 <- data[bootstrapIndices[[treatments[1]]], ]
        bootData2 <- data[bootstrapIndices[[treatments[2]]], ]
      }
      
      mi_matrix1 <- build.mim(bootData1, estimator = estimator)
      mi_matrix2 <- build.mim(bootData2, estimator = estimator)
      
      # Compare observed vs bootstrapped
      obs1_vs_boot2_higher <- obs1_vs_boot2_higher + (OBS_mi_matrix1 >= mi_matrix2)
      obs2_vs_boot1_higher <- obs2_vs_boot1_higher + (OBS_mi_matrix2 >= mi_matrix1)
      obs1_vs_boot2_lower <- obs1_vs_boot2_lower + (OBS_mi_matrix1 <= mi_matrix2)
      obs2_vs_boot1_lower <- obs2_vs_boot1_lower + (OBS_mi_matrix2 <= mi_matrix1)
      
      if(i %% 1000 == 0){
        dtime <- Sys.time() - StartTime
        rate <- i / as.numeric(dtime, units = "secs")
        remaining_iterations <- n_bootstrap - i
        eta_seconds <- remaining_iterations / rate
        if(eta_seconds > 3600) {
          eta_readable <- paste(round(eta_seconds/3600, 1), "hours")
        } else if(eta_seconds > 60) {
          eta_readable <- paste(round(eta_seconds/60, 1), "minutes") 
        } else {
          eta_readable <- paste(round(eta_seconds), "seconds")
        }
        cat("ETA:", eta_readable, "; (", round(i/n_bootstrap*100, 1), "%)\n")
      }
    }
  }else{
    if(workers == 0){
      workers <- detectCores() - 1
    }
    cl <- makeCluster(workers, outfile = "worker_output.txt")
    registerDoSNOW(cl)
    
    pb <- txtProgressBar(max = n_bootstrap, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    StartTime <- Sys.time()
    
    results <- foreach(i = 1:n_bootstrap, .combine = "+", 
                       .packages = c('minet'),
                       .export = c('StratifiedBootstrapIndices', 'TimepointCountData', 
                                   'OBS_mi_matrix1', 'OBS_mi_matrix2'),
                       .options.snow = opts) %dopar% {
                         gc()
                         set.seed(i)
                         if(!timeStrat){
                           bootIndex1 <- sample(1:nrow(group1_data), size = nrow(group1_data), replace = TRUE)
                           bootIndex2 <- sample(1:nrow(group2_data), size = nrow(group2_data), replace = TRUE)
                           
                           bootData1 <- group1_data[bootIndex1, ]
                           bootData2 <- group2_data[bootIndex2, ]
                         }else{
                           bootstrapIndices <- StratifiedBootstrapIndices(metadata, timeStratifiedMetadata, treatmentColumn)
                           bootData1 <- data[bootstrapIndices[[treatments[1]]], ]
                           bootData2 <- data[bootstrapIndices[[treatments[2]]], ]
                         }
                         
                         mi_matrix1 <- build.mim(bootData1, estimator = estimator)
                         mi_matrix2 <- build.mim(bootData2, estimator = estimator)
                         
                         # Return 3D array: [,,1] = obs1 >= boot2, [,,2] = obs2 >= boot1
                         array(c((OBS_mi_matrix1 >= mi_matrix2) * 1,
                                 (OBS_mi_matrix2 >= mi_matrix1) * 1,
                                 (OBS_mi_matrix1 <= mi_matrix2) * 1,
                                 (OBS_mi_matrix2 <= mi_matrix1) * 1), 
                               dim = c(nrow(mi_matrix1), ncol(mi_matrix1), 4))
                       }
    
    close(pb)
    obs1_vs_boot2_higher <- results[,,1]
    obs2_vs_boot1_higher <- results[,,2]
    
    obs1_vs_boot2_lower <- results[,,3]
    obs2_vs_boot1_lower <- results[,,4]
    
    stopCluster(cl)
    print(Sys.time()-StartTime)
  }
  
  # Calculate p-values
  p_matrix1_higher <- obs1_vs_boot2_higher / n_bootstrap
  p_matrix1_lower <- obs1_vs_boot2_lower / n_bootstrap
  p_matrix2_higher <- obs2_vs_boot1_higher / n_bootstrap
  p_matrix2_lower <- obs2_vs_boot1_lower / n_bootstrap
  
  p_matrix1_twotail <- pmin(pmin(p_matrix1_higher, p_matrix1_lower) * 2, 1)
  p_matrix2_twotail <- pmin(pmin(p_matrix2_higher, p_matrix2_lower) * 2, 1)
  
  dimnames(p_matrix1_twotail) <- dimnames(OBS_mi_matrix1)
  dimnames(p_matrix2_twotail) <- dimnames(OBS_mi_matrix2)
  
  # Create results dataframe
  ResultsDF1 <- p_matrix1_twotail %>%
    as.data.frame() %>%
    rownames_to_column("protein1") %>%
    pivot_longer(-protein1, names_to = "protein2", values_to = "pval_1vs2")
  
  ResultsDF2 <- p_matrix2_twotail %>%
    as.data.frame() %>%
    rownames_to_column("protein1") %>%
    pivot_longer(-protein1, names_to = "protein2", values_to = "pval_2vs1")
  
  OBS_MI1 <- OBS_mi_matrix1 %>%
    as.data.frame() %>%
    rownames_to_column("protein1") %>%
    pivot_longer(-protein1, names_to = "protein2", values_to = "OBS_MI1")
  
  OBS_MI2 <- OBS_mi_matrix2 %>%
    as.data.frame() %>%
    rownames_to_column("protein1") %>%
    pivot_longer(-protein1, names_to = "protein2", values_to = "OBS_MI2")
  
  # Combine results
  ResultsDF <- merge(ResultsDF1, ResultsDF2, by = c("protein1", "protein2"))
  ResultsDF$OBS_MI1 <- OBS_MI1$OBS_MI1
  ResultsDF$OBS_MI2 <- OBS_MI2$OBS_MI2
  ResultsDF$OBS_Dif <- ResultsDF$OBS_MI1 - ResultsDF$OBS_MI2
  
  # Add q-values
  ResultsDF$Qval_1vs2 <- qvalue(ResultsDF$pval_1vs2)$qvalues
  ResultsDF$Qval_2vs1 <- qvalue(ResultsDF$pval_2vs1)$qvalues
  
  # Filter to unique pairs
  ResultsDF <- ResultsDF[ResultsDF$protein1 != ResultsDF$protein2 & ResultsDF$protein1 < ResultsDF$protein2, ]
  
  gc()
  return(ResultsDF)
}

CategorizeResults <- function(results,
                              qThreshold = 0.05,
                              effectThreshold = 0.01) {
  
  results$category <- dplyr::case_when(
    results$Qval_1vs2 < qThreshold &
      results$Qval_2vs1 < qThreshold &
      abs(results$OBS_Dif) >= effectThreshold &
      results$OBS_Dif > 0 ~ "Significant Up",
    results$Qval_1vs2 < qThreshold &
      results$Qval_2vs1 < qThreshold &
      abs(results$OBS_Dif) >= effectThreshold &
      results$OBS_Dif < 0 ~ "Significant Down",
    TRUE ~ "Non-significant"
  )
  
  n_up <- sum(results$category == "Significant Up")
  n_down <- sum(results$category == "Significant Down")
  message(paste("Significant Up:", n_up, "| Significant Down:", n_down,
                "| Total pairs:", nrow(results)))
  
  return(results)
}