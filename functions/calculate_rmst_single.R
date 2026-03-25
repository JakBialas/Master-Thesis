calculate_rmst_single <- function(control_data, regimen_data, tau = 1460) {

  data_combined <- rbind(
    data.frame(time = control_data$time, event = control_data$event, arm = 0),
    data.frame(time = regimen_data$time, event = regimen_data$event, arm = 1)
  )
  

  result <- rmst2(time = data_combined$time,
                  status = data_combined$event,
                  arm = data_combined$arm,
                  tau = tau)
  

  return(result$unadjusted.result[1])
}
