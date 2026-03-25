w_KS_stat <- function(RCB_ctrl, RCB_exp, scale=0, ifplot=T, plots_folder=".", prj_title="Project") {
# Estimate Treatment Efficacy Score (TES) using two sample weighted Kolmogorov-Smirnov test
# by comparing RCB score from experimental treatment cohort vs control cohort and Area Beetwen CDF curves (ABC)
#
# Input:
# RCB_ctrl - vector of RCB values in control cohort
# RCB_exp - vector of RCB values in experimental treatment cohort
# scale - weight function parameter. Higher value gives higher weights to low RCB scores (non-negative value)
# ifplot - if plot should be generated and saved (logical)
# plots_folder - plot saving destination
# prj_title - name of the project to include in plot file name
#
# Output: TES and ABC(area between weighted eCDF curves)
# 
# Author: Michal Marczyk (michal.marczyk@yale.edu; michal.marczyk@polsl.edu)
  
  # create weighted eCDF
  ecdf_ctrl <- w_eCDF(RCB_ctrl, scale=scale)
  ecdf_exp <- w_eCDF(RCB_exp, scale=scale)
  ecdf_diff <- function(x){ecdf_exp(x) - ecdf_ctrl(x)}
  
  # calculate Treatment Efficacy score (TES) as an area between two eCDF curves
  x <- sort(unique(c(0,RCB_ctrl,RCB_exp)))
  ecdf_diff_calc <- ecdf_diff(x)
  TES <- sum(ecdf_diff_calc[-length(ecdf_diff_calc)] * (x[-1] - x[-length(x)]))/max(x)

  if (ifplot){
    require(gridExtra)
    require(ggplot2)
    
    data_plot <- data.frame(x)
    data_plot$ecdf_ctrl <- ecdf_ctrl(x)
    data_plot$ecdf_exp <- ecdf_exp(x)
    data_plot$ecdf_diff <- ecdf_diff_calc
    
    p1 <- ggplot(data_plot, aes(x=x, y=ecdf_ctrl)) + geom_step(col="red") + 
      geom_step(aes(y=ecdf_exp),col="blue") + theme_bw() + theme(plot.title=element_text(hjust=0.5)) +
      labs(x="RCB score", y="Weighted empirical CDF", title="Experimental (blue) vs control (red) treatment") + ylim(c(0,1))
    p2 <- ggplot(data_plot, aes(x=x, y=ecdf_diff)) + geom_step(col="green") + theme_bw() +
      labs(x="RCB score", y="Running Treatment Efficacy Score", title=paste0("TES=",signif(TES,2))) +
      theme(plot.title=element_text(hjust=0.5)) + geom_hline(yintercept=0, col="black")
    pdf(file=paste0(plots_folder,"/wKS_",prj_title,".pdf"), width=8, height=6)
    gridExtra::grid.arrange(p1,p2,nrow=2)
    dev.off()
  }
  
  return(TES)
}
