augment_data_random <- function(df, arm, subtype, percent = 25, threshold = 0.5) {
  dataf <- df[df$arm == arm & df$Subtype == subtype, ]
  
  if (nrow(dataf) == 0) {
    warning(paste("No data for the group", arm, "and subtype", subtype))
    return(data.frame())
  }
  
  data <- data.frame(
    time = dataf$efs.time,
    event = dataf$efs.ind,
    feature1 = dataf$rcbindex
  )
  
  log_mean_f1 <- mean(log(data$feature1[data$feature1 > 0]), na.rm = TRUE)
  log_sd_f1 <- sd(log(data$feature1[data$feature1 > 0]), na.rm = TRUE)
  
  log_mean_time <- mean(log(as.numeric(data$time[data$time > 0])), na.rm = TRUE)
  log_sd_time <- sd(log(as.numeric(data$time[data$time > 0])), na.rm = TRUE)
  
  gdata <- round(nrow(data) * (percent / 100))
  
  if (gdata > 0) {
    new_obs_lognorm <- data.frame(
      time = rlnorm(gdata, meanlog = log_mean_time, sdlog = log_sd_time),
      event = runif(gdata) < threshold,
      feature1 = rlnorm(gdata, meanlog = log_mean_f1, sdlog = log_sd_f1)
    )
    
    new_obs_lognorm$feature1 <- pmax(pmin(new_obs_lognorm$feature1, 5), 0)
    
    new_obs_lognorm$source <- "Generated"
    data$source <- "Original"
    
    combined_data <- rbind(data[, c("time", "event", "feature1", "source")],
                           new_obs_lognorm[, c("time", "event", "feature1", "source")])
    
    return(combined_data)
  } else {
    return(data)
  }
}
