generate_smoteCox_data <- function(df, arm, subtype, percent) {
  library(survival)
  library(smotefamily)
  
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
    warning("There are not enough classes in the 'event' variable to perform SMOTE.")
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
  
  cox_model <- coxph(Surv(time, event) ~ feature1, data = data)
  
  if (sum(data$event) == 0) {
    stop("No events in data, unable to calculate lambda 0")
  }
  lambda0 <- 1 / mean(data$time[data$event == 1])
  
  if (percent > 0) {
    num_synthetic <- round(nrow(data) * percent / 100)
    
    Y <- factor(data$event)
    X <- data.frame(feature1 = data$feature1, dummy = 0)
    smote_result <- SMOTE(X, data$event, K = K_value, dup_size = 1)
    
    synthetic_data <- smote_result$syn_data
    
    if (is.null(synthetic_data) || nrow(synthetic_data) == 0) {
      warning("No synthetic data. Returning original data.")
      data$source <- "Original"
      return(data)
    }
    
    synthetic_data$event <- as.numeric(as.character(synthetic_data$class))
    synthetic_data$class <- NULL
    
    risk_scores <- predict(cox_model, newdata = synthetic_data, type = "risk")
    synthetic_data$time <- rexp(nrow(synthetic_data), rate = lambda0 * risk_scores)
    
    synthetic_data$feature1 <- pmax(pmin(synthetic_data$feature1, 5), 0)
    
    synthetic_data$source <- "Generated"
    data$source <- "Original"
    
    combined_data <- rbind(
      data[, c("time", "feature1", "event", "source")],
      synthetic_data[, c("time", "feature1", "event", "source")]
    )
    
    return(head(combined_data, n = nrow(data) + num_synthetic))
    
  } else {
    data$source <- "Original"
    return(data)
  }
}
