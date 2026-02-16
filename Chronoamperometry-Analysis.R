# ==============================================================================
# Chronoamperometry Data Analysis 
# Automated extraction, statistical analysis, and power-law fitting
# ==============================================================================

library(readxl)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(openxlsx)
library(ggforce)

# ------------------------------------------------------------------------------
# 1. Data Import and Cleaning
# ------------------------------------------------------------------------------
# Place your data file in a folder named 'data' within your repository.
data_path <- "./data/PalmSens_Output.xlsx" 
raw_data <- read_excel(data_path)

# Standardize column names (assuming alternating Time and Current columns)
orig_names <- colnames(raw_data)
new_names <- character(length(orig_names))

for (i in seq_along(orig_names)) {
  if (i %% 2 == 1) {  
    # Odd columns: TIME
    new_names[i] <- paste0("TIME_", orig_names[i])
  } else {            
    # Even columns: CURRENT - inherit the identifier from the preceding Time column
    prev_name <- new_names[i-1]  
    current_name <- sub("TIME_", "CURRENT_", prev_name)  
    new_names[i] <- current_name
  }
}

colnames(raw_data) <- new_names

# Remove the first row if it contains units or metadata, and convert to numeric
raw_data <- raw_data[-1, ]

col_current <- seq(2, ncol(raw_data), by = 2)
raw_data[col_current] <- lapply(raw_data[col_current], as.numeric)

col_time <- seq(1, ncol(raw_data), by = 2)
raw_data[col_time] <- lapply(raw_data[col_time], as.numeric)

# ------------------------------------------------------------------------------
# 2. Extract Specific Data Points (Current at 5 seconds)
# ------------------------------------------------------------------------------
time_cols <- grep("TIME_", names(raw_data), value = TRUE)  
current_cols <- grep("CURRENT_", names(raw_data), value = TRUE)
results_5s <- data.frame()

for (i in seq_along(time_cols)) {
  time_vec <- raw_data[[time_cols[i]]]
  current_vec <- raw_data[[current_cols[i]]]
  
  # Find the index closest to 5 seconds
  idx_5s <- which.min(abs(time_vec - 5))
  current_5s <- current_vec[idx_5s]
  
  # Determine experimental condition based on column string
  condition <- if(grepl("Positivo", current_cols[i], ignore.case = TRUE)) {
    "Positive"
  } else if(grepl("Negativo", current_cols[i], ignore.case = TRUE)) {
    "Negative"
  } else {
    "Ultrapure Water"
  }
  
  results_5s <- rbind(results_5s, data.frame(
    Measurement = i,
    Time = time_vec[idx_5s],
    Current = current_5s,
    Condition = condition
  ))
}

results_5s$Condition <- as.factor(results_5s$Condition)

# ------------------------------------------------------------------------------
# 3. Statistical Analysis & Visualization (5s Point)
# ------------------------------------------------------------------------------
stats_summary <- results_5s %>%
  group_by(Condition) %>%
  summarise(
    Mean = mean(Current, na.rm = TRUE),
    Standard_Deviation = sd(Current, na.rm = TRUE)
  )

# Boxplot of Current at 5s
plot_5s <- ggplot(results_5s, aes(x = Condition, y = Current, fill = Condition)) +
  geom_boxplot(alpha = 0.8) +
  labs(
    title = "Current Distribution at t = 5 Seconds",
    x = "Experimental Condition",
    y = "Current (µA)"
  ) +
  scale_fill_manual(values = c("Positive" = "#E41A1C", "Negative" = "#377EB8", "Ultrapure Water" = "#4DAF4A")) +
  stat_compare_means(method = "t.test", label.y = max(results_5s$Current) * 1.1) + 
  theme_minimal()

print(plot_5s)

# ------------------------------------------------------------------------------
# 4. Power-Law Regression (Cottrell Equation & Anomalous Diffusion)
# Fitting model: I = a * t^b via logarithmic transformation: ln(I) = ln(a) + b*ln(t)
# ------------------------------------------------------------------------------
fit_log_linear <- function(time_val, current_val, sample_name) {
  
  # Filter data for the transient decay region (0.1 to 5 seconds)
  df_filtered <- data.frame(time = time_val, current = current_val) %>%
    filter(time > 0.1 & time <= 5 & !is.na(current) & !is.na(time) & current > 0)
  
  df_filtered <- df_filtered %>%
    mutate(log_time = log(time),
           log_current = log(current))
  
  # Linear model on log-transformed data
  model <- lm(log_current ~ log_time, data = df_filtered)
  
  log_a <- coef(model)[1]
  b <- coef(model)[2]
  a <- exp(log_a)
  r_squared <- summary(model)$r.squared
  
  # Generate fitted curve for visualization
  fitted_curve <- data.frame(time = seq(0.1, 5, length.out = 100)) %>%
    mutate(log_time = log(time),
           log_current = log_a + b * log_time,
           current = exp(log_current))
  
  p <- ggplot(df_filtered, aes(x = time, y = current)) +
    geom_point(color = "gray30", alpha = 0.5) +
    geom_line(data = fitted_curve, aes(x = time, y = current), color = "red", linewidth = 1) +
    labs(
      title = paste("Power-Law Fit (I = a*t^b):", sample_name),
      subtitle = bquote(R^2 == .(round(r_squared, 4)) ~ "|" ~ b == .(round(b, 4))),
      x = "Time (s)",
      y = "Current (µA)"
    ) +
    theme_minimal()
  
  return(list(params = data.frame(a = a, b = b, r2 = r_squared), plot = p))
}

# Execute fitting over all measurements
regression_results <- data.frame()
plot_list <- list()

for (i in seq_along(time_cols)) {
  t_data <- raw_data[[time_cols[i]]]
  i_data <- raw_data[[current_cols[i]]]
  
  cond <- results_5s$Condition[i] # Pulling condition mapped in step 2
  
  fit_output <- fit_log_linear(t_data, i_data, paste("Measurement", i, "-", cond))
  
  regression_results <- rbind(regression_results, 
                              data.frame(Condition = cond, 
                                         Measurement = i,
                                         a = fit_output$params$a, 
                                         b = fit_output$params$b,
                                         r2 = fit_output$params$r2))
  plot_list[[i]] <- fit_output$plot
}

# Boxplots for Regression Parameters
boxplot_a <- ggplot(regression_results, aes(x = Condition, y = a, fill = Condition)) +
  geom_boxplot(alpha = 0.8) +
  labs(title = "Pre-exponential Factor (a)", x = "Condition", y = "Parameter 'a'") +
  scale_fill_manual(values = c("Positive" = "#E41A1C", "Negative" = "#377EB8", "Ultrapure Water" = "#4DAF4A")) +
  theme_minimal()

boxplot_b <- ggplot(regression_results, aes(x = Condition, y = b, fill = Condition)) +
  geom_boxplot(alpha = 0.8) +
  labs(title = "Decay Exponent (b)", x = "Condition", y = "Parameter 'b'") +
  scale_fill_manual(values = c("Positive" = "#E41A1C", "Negative" = "#377EB8", "Ultrapure Water" = "#4DAF4A")) +
  theme_minimal()

combined_params <- ggarrange(boxplot_a, boxplot_b, ncol = 2, nrow = 1, common.legend = TRUE, legend = "bottom")
print(combined_params)

# ------------------------------------------------------------------------------
# 5. Advanced Visualizations
# ------------------------------------------------------------------------------
# Scatter plot clustering of exponents (b) vs. pre-exponential factors (a)
scatter_clustering <- ggplot(regression_results, aes(x = a, y = b, color = Condition)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_mark_ellipse(
    aes(fill = Condition, label = Condition), 
    inherit.aes = FALSE, x = regression_results$a, y = regression_results$b,
    show.legend = FALSE, alpha = 0.15
  ) +
  scale_color_manual(values = c("Positive" = "#E41A1C", "Negative" = "#377EB8", "Ultrapure Water" = "#4DAF4A")) +
  scale_fill_manual(values = c("Positive" = "#E41A1C", "Negative" = "#377EB8", "Ultrapure Water" = "#4DAF4A")) +
  labs(
    title = "Clustering of Diffusion Parameters",
    subtitle = "Decay exponent (b) vs Pre-exponential factor (a)",
    x = "Pre-exponential Factor (a)",
    y = "Decay Exponent (b)"
  ) +
  theme_minimal()

print(scatter_clustering)

# ------------------------------------------------------------------------------
# 6. Export Results
# ------------------------------------------------------------------------------
wb <- createWorkbook()

addWorksheet(wb, "5s_Results")
writeData(wb, sheet = "5s_Results", results_5s, rowNames = FALSE)

addWorksheet(wb, "Statistics")
writeData(wb, sheet = "Statistics", stats_summary, rowNames = FALSE)

addWorksheet(wb, "Curve_Analysis")
writeData(wb, sheet = "Curve_Analysis", regression_results, rowNames = FALSE)

export_path <- "./data/Consolidated_Results.xlsx"
saveWorkbook(wb, export_path, overwrite = TRUE)

message("Analysis complete. Results saved to: ", export_path)
