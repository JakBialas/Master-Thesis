RCB_ratio_fast <- function(RCB_ctrl, RCB_exp, sig=10^seq(-1,0.2,length.out=50), lam=10^seq(-3,0,length.out=50),
                      alpha=0.1, kernel_num=1000, n_pts=1000){
# Short version for calling within permutation test.
# 
# Author: Michal Marczyk (michal.marczyk@yale.edu; michal.marczyk@polsl.edu)
  
  require(densratio)
  
  # compute density ratio
  rat <- densratio(RCB_exp, RCB_ctrl, method="RuLSIF", verbose=F,
                   sigma=sig, lambda=lam, alpha=alpha, kernel_num=kernel_num)

  # Compute PDF and CDF
  RCBmax <- max(c(RCB_ctrl, RCB_exp))
  x_new <- seq(0,RCBmax,length.out=n_pts)
  w_hat <- rat$compute_density_ratio(x_new)
  cdf_hat <- cumsum(w_hat)
  cdf_hat <- cdf_hat/max(cdf_hat)

  # Calculate area under eCDF(Treatment Efficacy Score)
  AUC <- sum(cdf_hat[-length(cdf_hat)] * (x_new[-1] - x_new[-length(x_new)]))/RCBmax  #step function method
  AUC
}


