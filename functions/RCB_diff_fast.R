RCB_diff_fast <- function(RCB_ctrl, RCB_exp, sigma=1, lambda=1, RCBmax=5, n_pts=1000){
# Estimate Treatment Efficacy Score (TES) by comparing densities of RCB in experimental treatment cohort vs control cohort

  require(pcg)

  n1 <- length(RCB_exp)
  n2 <- length(RCB_ctrl)
  
  RCB_all <- c(RCB_exp, RCB_ctrl)
  RCBmax <- max(RCB_all)
  # rm(RCB_ctrl, RCB_exp)
  n <- n1+n2;
  kernel_num <- min(1000,n); # Number of kernel bases
  RCB_tmp <- seq(0,RCBmax,length.out=n_pts)
  
  center_index_tmp <- sample.int(n)
  center_index <- center_index_tmp[1:kernel_num]
  C <- RCB_all[center_index] 
  XXsum <- RCB_all^2
  XC_dist2 <- matrix(rep(XXsum,each=kernel_num), nrow=kernel_num,byrow=RCB_tmp) + 
    t(matrix(rep(XXsum[center_index],each=n), nrow=n, byrow=RCB_tmp)) - 2*as.matrix(C) %*% t(as.matrix(RCB_all))
  TC_dist2 <- matrix(rep(RCB_tmp^2,each=kernel_num), nrow=kernel_num,byrow=RCB_tmp) + 
    t(matrix(rep(XXsum[center_index],each=n_pts), nrow=n_pts, byrow=RCB_tmp)) - 2*as.matrix(C) %*% t(as.matrix(RCB_tmp))
  rm(XXsum,RCB_all)
  
  X1C_dist2 <- XC_dist2[,1:n1]
  X2C_dist2 <- XC_dist2[,(n1+1):n]
  CC_dist2 <- XC_dist2[,center_index]
  rm(XC_dist2)
  
  # Find density difference for optimal parameters
  H <- (sqrt(pi)*sigma)*exp(-CC_dist2/(4*sigma^2))
  h <- as.matrix(rowMeans(exp(-X1C_dist2/(2*sigma^2)))-rowMeans(exp(-X2C_dist2/(2*sigma^2))))
  thetah <- as.matrix(pcg(H+lambda*diag(kernel_num), h))

  # Compute Density difference and its CDF
  w_hat <- as.numeric(t(thetah) %*% exp(-TC_dist2/(2*sigma^2)))
  w_hat_pos <- w_hat
  w_hat_pos[w_hat<0] <- 0
  cdf_hat <- cumsum(w_hat_pos)
  cdf_hat <- cdf_hat/max(cdf_hat)
  
  # Calculate area under eCDF(Treatment Efficacy Score) with 95% confidence intervals
  TES <- sum(cdf_hat[-length(cdf_hat)] * (RCB_tmp[-1] - RCB_tmp[-length(RCB_tmp)]))/RCBmax  #step function method
 
  return(TES)
}