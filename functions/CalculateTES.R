CalculateTES <- function(Control, Regimen, regimen, subtype) {
 
  method <- "wKS" #"wKS","DensRatio" or "DensDiff"
  n_perm <- 10000
  n_cores <- 1
  scale <- -0.136 
  control_treatment <- Control
  experimental_treatment <- Regimen
  
  f_namen <- paste0(deparse(substitute(Control)), "_vs_", regimen,"_", subtype)
  # Definicja zmiennej f_name
  f_name <- paste0(deparse(substitute(Control)), "_vs_", deparse(substitute(Regimen)))
  
  if(!dir.exists(paste0("res/",f_name))) dir.create(paste0("res/",f_name), recursive=T)
  
  # Loading RCB score values (single column .txt files stored in data folder)
  RCB_ctrl <- as.numeric(control_treatment)
  RCB_exp <- as.numeric(experimental_treatment)
  cat("Control treatment cohort size:", length(RCB_ctrl), "\n")
  cat("Experimental treatment cohort size:", length(RCB_exp), "\n")
  
  # Calculating Treatment Efficacy Score and its p-value
  res_TES <- switch(method,
                    wKS = wKS(RCB_ctrl, RCB_exp, n_perm=n_perm, plots_folder=paste0("res/", f_name), prj_title=f_namen, n_cores=n_cores, scale=scale),
                    DensRatio = DensRatio(RCB_ctrl, RCB_exp, n_perm=n_perm, plots_folder=paste0("res/", f_name), prj_title=f_namen, n_cores=n_cores),
                    DensDiff = DensDiff(RCB_ctrl, RCB_exp, n_perm=n_perm, sig="auto", plots_folder=paste0("res/", f_name), prj_title=f_namen, n_cores=n_cores))
  
  return(res_TES)
}