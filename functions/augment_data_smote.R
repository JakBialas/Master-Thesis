generate_smote_data <- function(df, arm, subtype, percent) {
  dataf <- df[df$arm == arm & df$Subtype == subtype, ]
  
  if (nrow(dataf) == 0) {
    warning("No data after filtering")
    return(data.frame())
  }
  
  data <- data.frame(
    time = dataf$efs.time,
    event = dataf$efs.ind,
    feature1 = dataf$rcbindex
  )
  
  if (length(unique(data$event)) < 2) {
    warning("There are not enough classes in the 'event' variable to perform SMOTE")
    return(data.frame())
  }
  
  event_counts <- table(data$event)
  minority_class_count <- min(event_counts)
  
  if (minority_class_count < 3) {
    K_value <- 1
  } else if (minority_class_count < 4) {
    K_value <- 2
  } else {
    K_value <- 3
  }
  
  data$event <- as.factor(data$event)
  
  if (percent > 0) {
    num_synthetic <- round(nrow(data) * percent / 100)
    
    smote_result <- SMOTE(data[, c("feature1", "time")], 
                          data$event, 
                          K = K_value,              
                          dup_size = 1)
    
    synthetic_data <- smote_result$syn_data
    
    if (is.null(synthetic_data)) {
      warning("No synthetic data. Returning original data.")
      data$source <- "Original"
      return(data)
    }
    
    synthetic_data$event <- as.factor(synthetic_data$class)
    synthetic_data$class <- NULL
    synthetic_data$source <- "Generated"
    
    # Ograniczenia RCB
    synthetic_data$feature1 <- pmax(synthetic_data$feature1, 0)
    synthetic_data$feature1 <- pmin(synthetic_data$feature1, 5)
    
    data$source <- "Original"
    
    combined_data <- rbind(
      data[, c("time", "feature1", "event", "source")],
      synthetic_data[, c("time", "feature1", "event", "source")]
    )
    
    combined_data$event <- as.numeric(as.character(combined_data$event))
    combined_data$time <- as.numeric(as.character(combined_data$time))

    return(head(combined_data, n = nrow(data) + num_synthetic))
  } else {
    data$source <- "Original"
    data$event <- as.numeric(as.character(data$event))
    data$time <- as.numeric(as.character(data$time))
    return(data)
  }
}
