# ==============================================================================
# Multi-Condition Chronoamperometry Analysis (Dilution Series)
# Automated extraction, statistical ANOVA, and log-linear regression
# ==============================================================================

library(readxl)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(openxlsx)

# ------------------------------------------------------------------------------
# 1. Data Import and Cleaning
# ------------------------------------------------------------------------------
data_path <- "./data/PalmSens_Output.xlsx"
raw_data <- read_excel(data_path)

orig_names <- colnames(raw_data)
new_names <- character(length(orig_names))

for (i in seq_along(orig_names)) {
  if (i %% 2 == 1) {  
    new_names[i] <- paste0("TIME_", orig_names[i])
  } else {            
    prev_name <- new_names[i-1]  
    new_names[i] <- sub("TIME_", "CURRENT_", prev_name)  
  }
}

colnames(raw_data) <- new_names
raw_data <- raw_data[-1, ]

col_current <- seq(2, ncol(raw_data), by = 2)
raw_data[col_current] <- lapply(raw_data[col_current], as.numeric)

col_time <- seq(1, ncol(raw_data), by = 2)
raw_data[col_time] <- lapply(raw_data[col_time], as.numeric)

# ------------------------------------------------------------------------------
# 2. Data Extraction (Current at 5s)
# ------------------------------------------------------------------------------
time_cols <- grep("TIME_", names(raw_data), value = TRUE)  
current_cols <- grep("CURRENT_", names(raw_data), value = TRUE)
results_5s <- data.frame()

for (i in seq_along(time_cols)) {
  time_vec <- raw_data[[time_cols[i]]]
  current_vec <- raw_data[[current_cols[i]]]
  
  idx_5s <- which.min(abs(time_vec - 5))
  current_5s <- current_vec[idx_5s]
  
  # Map experimental conditions
  condition <- case_when(
    grepl("Basal", current_cols[i], ignore.case = TRUE) ~ "Basal",
    grepl("Controle -", current_cols[i], ignore.case = TRUE) ~ "Negative Control",
    grepl("1:2", current_cols[i]) ~ "1:2",
    grepl("1:8", current_cols[i]) ~ "1:8",
    grepl("1:32", current_cols[i]) ~ "1:32",
    grepl("1:128", current_cols[i]) ~ "1:128",
    grepl("1:1024", current_cols[i]) ~ "1:1024",
    TRUE ~ "Other"
  )
  
  results_5s <- rbind(results_5s, data.frame(
    Measurement = i,
    Time = time_vec[idx_5s],
    Current = current_5s,
    Condition = condition
  ))
}

# Define ordinal levels for correct plotting order
condition_order <- c("Negative Control", "1:2", "1:8", "1:32", "1:128", "1:1024", "Basal")
results_5s$Condition <- factor(results_5s$Condition, levels = condition_order)

# ------------------------------------------------------------------------------
# 3. Statistical Analysis & ANOVA Plot
# ------------------------------------------------------------------------------
stats_summary <- results_5s %>%
  group_by(Condition) %>%
  summarise(
    Mean = mean(Current, na.rm = TRUE),
    Std_Dev = sd(Current, na.rm = TRUE)
  )

anova_res <- aov(Current ~ Condition, data = results_5s)

plot_anova <- ggplot(results_5s, aes(x = Condition, y = Current, fill = Condition)) +
  geom_boxplot(alpha = 0.8) +
  labs(
    title = "Current Distribution at 5 Seconds",
    x = "Titration Condition",
    y = "Current (µA)"
  ) +
  scale_fill_manual(values = colorRampPalette(c("#E5F5E0", "#31A354"))(nlevels(results_5s$Condition))) + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(plot_anova)

# ------------------------------------------------------------------------------
# 4. Power-Law Regression (I = a*t^b)
# ------------------------------------------------------------------------------
fit_log_linear <- function(time_val, current_val) {
  df_filtered <- data.frame(time = time_val, current = current_val) %>%
    filter(time > 0, time <= 5, !is.na(current), !is.na(time), current > 0)
  
  df_filtered <- df_filtered %>%
    mutate(log_time = log(time), log_current = log(current))
  
  model <- lm(log_current ~ log_time, data = df_filtered)
  
  a <- exp(coef(model)[1])
  b <- coef(model)[2]
  r_squared <- summary(model)$r.squared
  
  return(list(a = a, b = b, r2 = r_squared))
}

regression_results <- data.frame()

for (i in seq_along(time_cols)) {
  if (!time_cols[i] %in% names(raw_data) || !current_cols[i] %in% names(raw_data)) next
  
  t_data <- raw_data[[time_cols[i]]]
  i_data <- raw_data[[current_cols[i]]]
  cond <- results_5s$Condition[i] # Inherit mapped condition
  
  # Defensive programming: tryCatch prevents loop from crashing on bad fits
  tryCatch({
    params <- fit_log_linear(t_data, i_data)
    regression_results <- rbind(regression_results,
                                data.frame(Condition = cond,
                                           a = params$a,
                                           b = params$b,
                                           r2 = params$r2,
                                           stringsAsFactors = FALSE))
  }, error = function(e) {
    message("Skipping fit for condition ", cond, " due to error: ", e$message)
  })
}

regression_results$Condition <- factor(regression_results$Condition, levels = condition_order)

# Regression Parameter Plots
plot_a <- ggplot(regression_results, aes(x = Condition, y = a, fill = Condition)) +
  geom_boxplot() +
  labs(title = "Pre-exponential Factor (a)", x = "Condition", y = "Value (a)") +
  scale_fill_manual(values = colorRampPalette(c("#FEE0D2", "#DE2D26"))(nlevels(regression_results$Condition))) +
  theme_minimal() + theme(axis.text.x = element_blank())

plot_b <- ggplot(regression_results, aes(x = Condition, y = b, fill = Condition)) +
  geom_boxplot() +
  labs(title = "Decay Exponent (b)", x = "Condition", y = "Value (b)") +
  scale_fill_manual(values = colorRampPalette(c("#FEE0D2", "#DE2D26"))(nlevels(regression_results$Condition))) +
  theme_minimal() + theme(axis.text.x = element_blank())

ggarrange(plot_a, plot_b, ncol = 2, common.legend = TRUE, legend = "bottom")

# ------------------------------------------------------------------------------
# 5. Export Results
# ------------------------------------------------------------------------------
wb <- createWorkbook()
addWorksheet(wb, "5s_Results")
writeData(wb, sheet = "5s_Results", results_5s, rowNames = FALSE)

addWorksheet(wb, "Statistics")
writeData(wb, sheet = "Statistics", stats_summary, rowNames = FALSE)

addWorksheet(wb, "Curve_Analysis")
writeData(wb, sheet = "Curve_Analysis", regression_results, rowNames = FALSE)

export_path <- "./data/Titration_Results.xlsx"
saveWorkbook(wb, export_path, overwrite = TRUE)
message("Analysis complete. Results saved to: ", export_path)