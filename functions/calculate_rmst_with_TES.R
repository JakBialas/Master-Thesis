calculate_rmst_with_TES <- function(df, TES_data,
                                    regimens_HER2 = c("Regimen6", "Regimen8", "Regimen9", "Regimen10", "Regimen11"),
                                    regimens_TNBC = c("Regimen1", "Regimen2", "Regimen3", "Regimen4", "Regimen5", "Regimen6", "Regimen7")) {
  df$efs.time <- as.numeric(df$efs.time)
  df <- na.omit(df)

  subtypes <- c("HER2+", "TNBC")

  results <- data.frame(Subtype = character(), Regimen = character(),
                        Arm0_RMST = numeric(), Arm1_RMST = numeric(), Delta_RMST = numeric())

  for (subtype in subtypes) {
    df_filtered <- df %>% filter(Subtype == subtype)
    regimens <- if (subtype == "HER2+") regimens_HER2 else regimens_TNBC

    for (regimen in regimens) {
      dataf_control <- df_filtered %>% filter(arm == "Control")
      dataf_regimen <- df_filtered %>% filter(arm == regimen)
      dataf <- rbind(dataf_control, dataf_regimen)
      dataf$arm <- ifelse(dataf$arm == "Control", 0, 1)

      data <- data.frame(
        time = dataf$efs.time,
        event = dataf$efs.ind,
        arm = dataf$arm
      )

      tau <- 1460

      r <- rmst2(data$time, data$event, data$arm, tau = tau)

      results <- rbind(results, data.frame(
        Subtype = subtype,
        Regimen = regimen,
        Arm0_RMST = r$RMST.arm0$rmst[1],
        Arm1_RMST = r$RMST.arm1$rmst[1],
        Delta_RMST = r$unadjusted.result[1]
      ))
    }
  }

  results <- results %>% left_join(TES_data, by = c("Regimen", "Subtype"))
  return(results)
}