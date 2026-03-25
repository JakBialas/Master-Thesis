generate_smote_data_2Ev <- function(df, arm, subtype, percent) {
  dataf <- df[df$arm == arm & df$Subtype == subtype, ]
  
  if (nrow(dataf) == 0) {
    warning("No data after filtering")
    return(data.frame())
  }
  
  if (any(is.na(dataf$rcbindex)) || any(is.na(dataf$efs.time))) {
    stop("Missing values in input data")
  }
  
  base_data <- data.frame(
    time = dataf$efs.time,
    event = dataf$efs.ind,
    feature1 = dataf$rcbindex,
    source = "Original"
  )

  if (percent <= 0) {
    return(base_data)
  }
  all_data <- base_data
  
  for (event_val in c(0, 1)) {
    data_event <- base_data[base_data$event == event_val, ]
    
    K_val <- if (nrow(data_event) < 3) 1 else if (nrow(data_event) < 4) 2 else 3
    
    smote_res <- tryCatch({
      SMOTE(data_event[, c("feature1", "time")],
            as.factor(rep(1, nrow(data_event))),
            K = K_val, dup_size = 1)
    }, error = function(e) {
      warning(paste("Error in SMOTE for event =", event_val))
      return(NULL)
    })
    
    if (is.null(smote_res) || is.null(smote_res$syn_data)) next
    
    syn <- smote_res$syn_data
    syn$event <- event_val
    syn$source <- "Generated"
    
    syn$feature1 <- pmax(pmin(syn$feature1, 5), 0)
    syn$time <- syn$time
    
    synthetic_count <- round(nrow(data_event) * percent / 100)
    if (synthetic_count < nrow(syn)) {
      syn <- syn[1:synthetic_count, ]
    }
    
    syn <- syn[, c("time", "feature1", "event", "source")]
    all_data <- rbind(all_data, syn)
  }
  
  return(all_data)
}
