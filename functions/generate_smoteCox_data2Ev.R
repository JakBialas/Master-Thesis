generate_smoteCox_data2Ev <- function(df, arm, subtype, percent) {
  library(survival)
  library(smotefamily)
  
  dataf <- df[df$arm == arm & df$Subtype == subtype, ]
  
  if (nrow(dataf) == 0) {
    warning("Brak danych po filtracji")
    return(data.frame())
  }
  
  if (any(is.na(dataf$rcbindex)) || any(is.na(dataf$efs.time))) {
    stop("No value in 'rcbindex' or 'efs.time'")
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
  
  if (length(unique(base_data$event)) > 1 && sum(base_data$event) > 0) {
    cox_model <- coxph(Surv(time, event) ~ feature1, data = base_data)
    lambda0 <- 1 / mean(base_data$time[base_data$event == 1])
  } else {
    cox_model <- NULL
    lambda0 <- NULL
  }
  
  for (event_val in unique(base_data$event)) {
    data_event <- base_data[base_data$event == event_val, ]
    
    if (nrow(data_event) < 2) next
    
    # Ustalenie K dla SMOTE
    K_val <- if (nrow(data_event) < 3) 1 else if (nrow(data_event) < 4) 2 else 3
    
    smote_res <- tryCatch({
      # Dodajemy dummy do SMOTE, jesli mamy tylko jedna ceche
      X <- data.frame(feature1 = data_event$feature1, dummy = 0)
      Y <- factor(rep(1, nrow(data_event)))
      SMOTE(X, Y, K = K_val, dup_size = 1)
    }, error = function(e) {
      warning(paste("Blad SMOTE dla event =", event_val))
      return(NULL)
    })
    
    if (is.null(smote_res) || nrow(smote_res$syn_data) == 0) next
    
    syn <- smote_res$syn_data
    syn$feature1 <- pmax(pmin(syn$feature1, 5), 0)
    syn$event <- event_val
    syn$source <- "Generated"
    
    if (!is.null(cox_model) && !is.null(lambda0)) {
      risk_scores <- predict(cox_model, newdata = syn, type = "risk")
      syn$time <- rexp(nrow(syn), rate = lambda0 * risk_scores)
    } else {
      syn$time <- syn$time
    }
    
    synthetic_count <- round(nrow(data_event) * percent / 100)
    if (synthetic_count < nrow(syn)) syn <- syn[1:synthetic_count, ]
    
    syn <- syn[, c("time", "feature1", "event", "source")]
    all_data <- rbind(all_data, syn)
  }
  
  return(all_data)
}
