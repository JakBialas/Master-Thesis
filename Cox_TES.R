setwd("/.../.../...")
options(repos = c(CRAN = "https://cloud.r-project.org"))
tmp <- sapply(paste0("Rscripts/", list.files(path = "Rscripts")), source)

libs <- c("ggplot2", "densratio", "parallel", "reshape", "pcg", "gridExtra", "dplyr", "writexl", "tidyr", "smotefamily", "survRM2")
check_libs(libs)
library(smotefamily)
library(survRM2)

df <- read.csv("/.../.../.../BC_Data.csv")
df <- na.omit(df)

subtypes <- unique(df$Subtype)
regimens_list <- df %>%
  filter(arm != "Control") %>%
  group_by(Subtype) %>%
  summarise(Regimens = list(unique(arm)), .groups = "drop")

combs <- do.call(rbind, lapply(subtypes, function(subtype) {
  regimens_for_subtype <- regimens_list$Regimens[regimens_list$Subtype == subtype][[1]]
  expand.grid(Subtype = subtype, Percent = seq(0, 100, 1), Regimen = regimens_for_subtype, stringsAsFactors = FALSE)
}))

numCores <- 40

result_list <- mclapply(1:nrow(combs), function(i) {
  row <- combs[i, ]
  subtype <- row$Subtype
  percent <- row$Percent
  regimen <- row$Regimen
  
  cat("Processing: Subtype =", subtype, ", Percent =", percent, "%, Regimen =", regimen, "\n")
  flush.console()
  
  TES_values <- numeric(10)
  p_values <- numeric(10)
  rmst_values <- numeric(10)  

  for (j in 1:10) {
    cat("Running repetition", j, "for augmentation at:", percent, "%\n")
    flush.console()
    
    ControlData <- generate_cox_data(df, 'Control', subtype, percent)
    if (nrow(ControlData) == 0) {
      warning(paste("No data for the control group at", percent, "% augmentation"))
      next
    }
    Control <- ControlData$feature1

    RegimenData <- generate_cox_data(df, regimen, subtype, percent)
    if (nrow(RegimenData) == 0) {
      warning(paste("No data for the regiment", regimen, "at", percent, "% augmentation"))
      next
    }
    Regimen <- RegimenData$feature1

    res <- CalculateTES(Control, Regimen, regimen, subtype, percent)
    if (!is.null(res) && length(res) == 2) {
      TES_values[j] <- res[1]  # TES
      p_values[j] <- res[2]    # p-value
    } else {
      warning(paste("Error while calculating TES for regimen", regimen, "at", percent, "% augmentation"))
    }
    rmst_val <- calculate_rmst_single(ControlData, RegimenData)
    rmst_values[j] <- rmst_val
  }
  
  TES_avg <- mean(TES_values, na.rm = TRUE)
  TES_var <- var(TES_values, na.rm = TRUE)
  p_value_avg <- mean(p_values, na.rm = TRUE)
  Delta_RMST_avg <- mean(rmst_values, na.rm = TRUE)

  data.frame(
    Subtype = subtype,
    Regimen = regimen,
    Percent = percent,
    TES_avg = TES_avg,
    TES_var = TES_var,
    p_value_avg = p_value_avg,
    Delta_RMST = Delta_RMST_avg,
    stringsAsFactors = FALSE
  )
}, mc.cores = numCores)

results <- do.call(rbind, result_list)
print(results)
excelname <- "COXv2.xlsx"
write_xlsx(results, excelname)

if (!dir.exists("plots_COXv2")) dir.create("plots_COXv2")

for (subtype in unique(results$Subtype)) {
  for (regimen in unique(results$Regimen)) {
    subset_results <- results[results$Subtype == subtype & results$Regimen == regimen, ]
    if (nrow(subset_results) > 0) {
      
      TES_original <- subset_results$TES_avg[subset_results$Percent == 0]
      TES = subset_results$TES_avg
      TES_avg <- subset_results$TES_avg
      TES_var <- subset_results$TES_var
      p_value_avg <- subset_results$p_value_avg
      
      subset_results <- subset_results %>%
        mutate(TES_diff = TES - TES_original,
               TES_avg = TES_avg,
               TES_var = TES_var)

      lm_model_TES <- lm(TES_diff ~ Percent, data = subset_results)
      slope_TES <- coef(lm_model_TES)[2]
      intercept_TES <- coef(lm_model_TES)[1]

      slope_TES_formatted <- sprintf("%.10f", slope_TES)
      intercept_TES_formatted <- sprintf("%.10f", intercept_TES)

      lm_model_variance <- lm(TES_var ~ Percent, data = subset_results)
      slope_variance <- coef(lm_model_variance)[2]
      intercept_variance <- coef(lm_model_variance)[1]

      slope_variance_formatted <- sprintf("%.10f", slope_variance)
      intercept_variance_formatted <- sprintf("%.10f", intercept_variance)

      plot_TES <- ggplot(subset_results, aes(x = Percent)) +
        geom_line(aes(y = TES_diff), color = "blue", linewidth = 1) +
        geom_point(aes(y = TES_diff), color = "red", linewidth = 3) +
        geom_hline(yintercept = 0, color = "green", linewidth = 1, linetype = "dotted") +
        geom_abline(slope = slope_TES, intercept = intercept_TES, color = "purple", linetype = "dashed", size = 1) +
        labs(title = paste("(COX) TES Dependency, (", subtype, " ", regimen, ")"),
             subtitle = paste("Slope =", slope_TES_formatted, ", Intercept =", intercept_TES_formatted),
             x = "Augmentation percentage", y = "TES") +
        theme_light() +
        ylim(-0.2, 0.4)

      plot_variance <- ggplot(subset_results, aes(x = Percent)) +
        geom_line(aes(y = TES_var), color = "purple", linewidth = 1) +
        geom_point(aes(y = TES_var), color = "orange", linewidth = 3) +
        geom_abline(slope = slope_variance, intercept = intercept_variance, color = "darkblue", linetype = "dashed", size = 1) +
        labs(title = paste("(COX) TES Variance, (", subtype, " ", regimen, ")"),
             subtitle = paste("Slope =", slope_variance_formatted, ", Intercept =", intercept_variance_formatted),
             x = "Augmentation percentage", y = "TES Variance") +
        theme_light()

      plot_TES_filename <- paste0("plots_COXv2/TES_", subtype, "_", regimen, ".png")
      plot_variance_filename <- paste0("plots_COXv2/Variance_", subtype, "_", regimen, ".png")
      
      ggsave(plot_TES_filename, plot = plot_TES, width = 8, height = 6, dpi = 300)
      ggsave(plot_variance_filename, plot = plot_variance, width = 8, height = 6, dpi = 300)
    }
  }
}

######
augmentation_percents <- seq(0, 100, by = 5)

if (!dir.exists("plots_COXv2")) dir.create("plots_COXv2")

for (augmentation_percent in augmentation_percents) {
  results_subset <- results %>% filter(Percent == augmentation_percent)
  
  if (nrow(results_subset) == 0) {
    warning(paste("No data for", augmentation_percent, "% augmentation"))
    next
  }
  
  pearson_corr <- cor(results_subset$TES_avg, results_subset$Delta_RMST, method = "pearson")
  spearman_corr <- cor(results_subset$TES_avg, results_subset$Delta_RMST, method = "spearman")
  
  lm_model <- lm(Delta_RMST ~ TES_avg, data = results_subset)
  slope <- coef(lm_model)[2]
  intercept <- coef(lm_model)[1]
  
  plot <- ggplot(results_subset, aes(x = TES_avg, y = Delta_RMST, color = Subtype, label = Regimen)) +
    geom_point(size = 3) +
    geom_text(vjust = -0.5, hjust = 0.5, size = 3) +
    geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
    labs(title = paste("Delta RMST vs TES for", augmentation_percent, "% Augmentation (COX)"),
         x = "TES", y = "Delta RMST") +
    theme_light() +
    theme(legend.position = "top") +
    annotate("text",
             x = -Inf, y = Inf,
             hjust = -0.05, vjust = 1.2,
             label = paste0("Pearson R = ", round(pearson_corr, 2),
                            "\nSpearman R = ", round(spearman_corr, 2),
                            "\nSlope = ", round(slope, 2),
                            "\nIntercept = ", round(intercept, 2)),
             size = 4, color = "black")
  
  filename <- paste0("plots_COXv2/Delta_RMST_vs_TES_", augmentation_percent, "_percent.png")
  ggsave(filename, plot = plot, width = 8, height = 6, dpi = 300)
}
cat("The graph was saved in the 'plots' folder.\n")


#####
augment_levels <- unique(results$Percent)
correlation_results <- data.frame()

for (percent in augment_levels) {
  rmst_data <- results %>% filter(Percent == percent)
  
  pearson_r <- cor(rmst_data$TES_avg, rmst_data$Delta_RMST, method = "pearson")
  spearman_r <- cor(rmst_data$TES_avg, rmst_data$Delta_RMST, method = "spearman")
  
  lm_model <- lm(Delta_RMST ~ TES_avg, data = rmst_data)
  slope <- coef(lm_model)[2]
  intercept <- coef(lm_model)[1]
  
  correlation_results <- rbind(correlation_results, data.frame(
    Percent = percent,
    Pearson_R = pearson_r,
    Spearman_R = spearman_r,
    Slope = slope,
    Intercept = intercept
  ))
}

if (!dir.exists("plots_COXv2")) dir.create("plots_COXv2")

plot_metric <- function(data, y_col, y_label) {
  p <- ggplot(data, aes(x = Percent, y = .data[[y_col]])) +
    geom_line(color = "blue", size = 1.2) +
    geom_point(color = "red", size = 3) +
    geom_smooth(method = "lm", color = "black", linetype = "dashed", se = TRUE) +
    labs(title = paste(y_label, "vs Augmentation (COX)"),
         x = "Augmentation Percentage", y = y_label) +
    theme_light()
  return(p)
}

plot_pearson <- plot_metric(correlation_results, "Pearson_R", "Pearson R")
plot_spearman <- plot_metric(correlation_results, "Spearman_R", "Spearman R")
plot_slope <- plot_metric(correlation_results, "Slope", "Slope")
plot_intercept <- plot_metric(correlation_results, "Intercept", "Intercept")

library(gridExtra)
grid_plot <- grid.arrange(plot_pearson, plot_spearman, plot_slope, plot_intercept, ncol = 2)

ggsave("plots_COXv2/Correlation_and_Regression_vs_Augmentation_2x2.png", plot = grid_plot, width = 12, height = 8, dpi = 300)
######


cat("The plots were saved in the 'plots_COXv2' folder.\n")
flush.console()

summary_results <- results %>%
  group_by(Subtype, Regimen) %>%
  mutate(
    TES_original = TES_avg[Percent == 0],
    TES_diff = TES_avg - TES_original,
    TES_avg = TES_avg,
    TES_var = TES_var,
    p_value_avg = p_value_avg
  )

cat("\nSummary of results:\n")
print(summary_results %>% 
        select(Subtype, Regimen, Percent, TES_avg, TES_diff, TES_var, p_value_avg) %>% 
        arrange(Subtype, Regimen, Percent))

summary_with_corr <- summary_results %>%
  left_join(correlation_results, by = "Percent")

excelname <- "COXv2.xlsx"
write_xlsx(summary_with_corr, excelname)