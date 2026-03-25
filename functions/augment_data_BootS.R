augment_data_BootS <- function(df, arm, subtype, percent) {
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
  
  num_samples <- round(nrow(data) * percent / 100)
  
  if (num_samples == 0) {
    data$source <- "Original"
    return(data[, c("time", "feature1", "event", "source")])
  }
  
  bootstrap_indices <- sample(1:nrow(data), size = num_samples, replace = TRUE)
  
  new_samples <- data[bootstrap_indices, ]
  
  new_samples$source <- "Generated"
  data$source <- "Original"
  
  combined_data <- rbind(data[, c("time", "feature1", "event", "source")],
                         new_samples[, c("time", "feature1", "event", "source")])
  
  return(combined_data)
}
