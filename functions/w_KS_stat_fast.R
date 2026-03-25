w_KS_stat_fast <- function(RCB_ctrl,RCB_exp, scale=0) {
  
  # create weighted eCDFs and calculated their difference
  ecdf_ctrl <- w_eCDF(RCB_ctrl, scale=scale)
  ecdf_exp <- w_eCDF(RCB_exp, scale=scale)
  ecdf_diff <- function(x){ecdf_exp(x) - ecdf_ctrl(x)}
  
  # calculate Treatment Efficacy score (TES) as an area between two eCDF curves
  x <- sort(unique(c(0,RCB_ctrl,RCB_exp)))
  ecdf_diff_calc <- ecdf_diff(x)
  TES <- sum(ecdf_diff_calc[-length(ecdf_diff_calc)] * (x[-1] - x[-length(x)]))/max(x)
  
  return(TES)
}
