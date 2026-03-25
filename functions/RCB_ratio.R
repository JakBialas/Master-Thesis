RCB_ratio <- function(RCB_ctrl, RCB_exp, sig=10^seq(-1,0.2,length.out=50), lam=10^seq(-3,0,length.out=50), n_pts=10000,
                      alpha=0.1, kernel_num=1000, ifplot=F, plots_folder=".", prj_title="Project"){
# Estimate Treatment Efficacy Score (TES) by comparing densities of RCB in experimental treatment cohort vs control cohort
# Step1: estimate density ratio using kernel functions
# Step2: calculate eCDF of density ratio
# Step3: calculate TES as area under eCDF
#
# Input:
# RCB_ctrl - vector of RCB values in control cohort
# RCB_exp - vector of RCB values in experimental treatment cohort
# sig - kernel bandwidth (vector with numeric values > 0)
# lam - regularization parameter (vector with numeric values > 0)
# alpha - relative parameter for RuLSIF density estimation method (numeric from 0 to 1)
# kernel_num - no. of kernels (integer > 0)
# ifplot - if plot should be generated and saved (logical)
# plots_folder - plot saving destination
# prj_title - name of the project to include in plot file name
#
# Output: TES with 95% confidence intervals
# 
# Author: Michal Marczyk (michal.marczyk@yale.edu; michal.marczyk@polsl.edu)
  
  require(densratio)
  
  if (sig == "auto"){
    sig <- bw.nrd0(c(RCB_ctrl, RCB_exp)) 
  }
  
  RCBmax <- max(c(RCB_ctrl, RCB_exp),na.rm=T)
  n_pts <- 10000

  # compute density ratio
  rat <- densratio(RCB_exp, RCB_ctrl, method="RuLSIF", verbose=F,
                   sigma=sig, lambda=lam, alpha=alpha, kernel_num=kernel_num)
  if(ifplot) cat("Estimated parameters: lambda =",rat$lambda, "; bandwidth =", rat$kernel_info$sigma,"\n")
  
  # Compute PDF and CDF
  x_new <- seq(0,RCBmax,length.out=n_pts)
  w_hat <- rat$compute_density_ratio(x_new)
  cdf_hat <- cumsum(w_hat)
  cdf_hat <- cdf_hat/max(cdf_hat)
  
  # Calculate area under eCDF(Treatment Efficacy Score) with 95% confidence intervals
  AUC <- sum(cdf_hat[-length(cdf_hat)] * (x_new[-1] - x_new[-length(x_new)]))/RCBmax  #step function method

  # plot results
  if(ifplot){
    require(ggplot2)
    require(gridExtra)
    
    # plot density in individual groups
    data_plot <- data.frame(c(RCB_ctrl, RCB_exp))
    colnames(data_plot) <- "RCBscore"
    data_plot$Group <- c(rep("Control",length(RCB_ctrl)), rep("Experimental",length(RCB_exp)))
    p1 <- ggplot(data_plot,aes(col=Group, x=RCBscore)) + geom_density() + 
      theme_bw() + theme(legend.position="bottom") + labs(y="Density")

    # plot density ratio
    data_plot <- data.frame(cbind(x_new,w_hat, cdf_hat))
    p2 <- ggplot(data_plot,aes(x=x_new, y=w_hat)) + geom_point() + 
      labs(x="RCBscore", y="Density ratio\n(Experimental vs Control)", title=paste0("TES =",signif(AUC,2))) + 
      geom_hline(yintercept=1) + theme_bw() + theme(plot.title = element_text(hjust=0.5))
    
    # save figure
    pdf(file=paste0(plots_folder,"/DensRatio_",prj_title,".pdf"), width=8, height=6)
    gridExtra::grid.arrange(p1,p2,nrow=2)
    dev.off()
  }
  res <- list()
  res[["TES"]] <- AUC   #c(AUC,AUC_CI)
  # names(res[["TES"]]) <- c("TES","LOW_CI", "HIGH_CI")
  res[["Density_ratio"]] <- rat$compute_density_ratio
  res[["Parameters"]] <- c(rat$alpha, rat$lambda, rat$kernel_info$sigma, rat$kernel_info$kernel_num)
  names(res[["Parameters"]]) <- c("Alpha","Lambda","Sigma","Kernel_num")
  res
}


