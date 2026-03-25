generate_cox_data <- function(df, arm, subtype, percent = 0) {
  library(dplyr)
  library(survival)
  
  req_cols <- c("arm", "Subtype", "efs.time", "efs.ind", "rcbindex")
  stopifnot(all(req_cols %in% names(df)))
  
  data <- df %>%
    filter(arm == !!arm, Subtype == !!subtype) %>%
    transmute(time     = efs.time,
              event    = efs.ind,
              feature1 = rcbindex)
  
  if (nrow(data) == 0) {
    warning("No data to filter")
    return(tibble::tibble())
  }
  if (length(unique(data$event)) < 2) {
    warning("Not enough classes for modeling")
    return(tibble::tibble())
  }
  
  cox_model <- coxph(Surv(time, event) ~ feature1, data = data)
  if (sum(data$event) == 0)
    stop("There are no events (event == 1); ?0 cannot be calculated")
  lambda0 <- 1 / mean(data$time[data$event == 1])
  
  data <- mutate(data, source = "Original")
  if (percent <= 0) return(data)
  
  n_syn <- round(nrow(data) * percent / 100)
  if (n_syn == 0) return(data)
  

  log_mean_f1 <- mean(log(data$feature1[data$feature1 > 0]), na.rm = TRUE)
  log_sd_f1   <- sd(log(data$feature1[data$feature1 > 0]), na.rm = TRUE)
  
  synthetic_feature1 <- rlnorm(n_syn, meanlog = log_mean_f1, sdlog = log_sd_f1)
  synthetic_feature1 <- pmin(pmax(synthetic_feature1, 0), 5)

  #synthetic_feature1 <- rnorm(n_syn, mean = mean(data$feature1), sd = sd(data$feature1)) #rnorm(n_syn, mean = 2.5, sd = 1.5) #rnorm(n_syn, mean = mean(data$feature1), sd = sd(data$feature1)) #runif(n_syn, min = 0, max = 5)
  #synthetic_feature1 <- pmin(pmax(synthetic_feature1, 0), 5)
  
  raw_risk <- predict(cox_model,
                      newdata = tibble(feature1 = synthetic_feature1),
                      type = "risk")
  synthetic_time  <- rexp(n_syn, rate = lambda0 * raw_risk)
  #synthetic_event <- as.integer(runif(n_syn) < 0.5) #  synthetic_event <- as.integer(runif(n_syn) < 0.5) #  synthetic_event <- rep(1L, n_syn)
   
  prop_event <- mean(data$event == 1)
  synthetic_event <- as.integer(runif(n_syn) < prop_event)

  synthetic <- tibble::tibble(
    time     = synthetic_time,
    event    = synthetic_event,
    feature1 = synthetic_feature1,
    source   = "Generated"
  )
  
  bind_rows(data, synthetic)
  combined_data <- bind_rows(data, synthetic)
  return(combined_data)
}
