augment_data_BootSwLN <- function(df, arm, subtype, percent, noise_sd_factor = 0.1) {
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
  
  log_mean_feature1 <- mean(log(data$feature1[data$feature1 > 0]), na.rm = TRUE)
  log_sd_feature1 <- sd(log(data$feature1[data$feature1 > 0]), na.rm = TRUE) * noise_sd_factor
  
  log_mean_time <- mean(log(data$time[data$time > 0]), na.rm = TRUE)
  log_sd_time <- sd(log(data$time[data$time > 0]), na.rm = TRUE) * noise_sd_factor
  
  noise_feature1 <- rlnorm(nrow(new_samples), meanlog = 0, sdlog = log_sd_feature1)
  noise_time <- rlnorm(nrow(new_samples), meanlog = 0, sdlog = log_sd_time)
  
  new_samples$feature1 <- new_samples$feature1 * noise_feature1
  new_samples$time <- new_samples$time * noise_time
  
  new_samples$feature1 <- pmax(pmin(new_samples$feature1, 5), 0)
  
  new_samples$source <- "Generated"
  data$source <- "Original"
  
  combined_data <- rbind(data[, c("time", "feature1", "event", "source")],
                         new_samples[, c("time", "feature1", "event", "source")])
  
  return(combined_data)
}
