# Parametric T-test for Event Study Analysis
# This script performs parametric t-tests to analyze abnormal returns and
# cumulative abnormal returns (CAR) around event dates.

# -----------------------------------------------------------------------------
# Dependencies
# -----------------------------------------------------------------------------
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  dplyr,      # Data manipulation
  writexl     # Excel output
)

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
# Set default options
options(stringsAsFactors = FALSE)
options(scipen = 999)

# Define constants
EVENT_WINDOWS <- list(
  day_0    = c(0, 0),
  window_1 = c(-1, 1),
  window_2 = c(-2, 2)
)

# -----------------------------------------------------------------------------
# Helper Functions
# -----------------------------------------------------------------------------
#' Create Directory if Not Exists
#' @param path Character string specifying directory path
#' @return NULL
create_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
}

#' Calculate Abnormal Returns for a Single Event
#' @param data Dataframe containing event data
#' @return Dataframe with calculated abnormal returns
calculate_abnormal_returns <- function(data) {
  # Split into estimation and event windows
  estimation_data <- data %>% filter(Relative_trading_day < 0)
  event_window    <- data %>% filter(Relative_trading_day >= 0)
  
  # Check for sufficient data
  if (nrow(estimation_data) == 0 || nrow(event_window) == 0) {
    warning("Insufficient data for event ", unique(data$Event_ID))
    return(data)
  }
  
  tryCatch({
    # Fit market model
    market_model <- lm(Stock_return ~ Market_return + World_return + FX_weighted_pct_change, 
                      data = estimation_data)
    
    # Calculate abnormal returns
    estimation_data <- estimation_data %>%
      mutate(
        predicted_return = predict(market_model, newdata = estimation_data),
        abnormal_return  = Stock_return - predicted_return
      )
    
    event_window <- event_window %>%
      mutate(
        predicted_return = predict(market_model, newdata = event_window),
        abnormal_return  = Stock_return - predicted_return
      )
    
    bind_rows(estimation_data, event_window)
    
  }, error = function(e) {
    warning("Error in market model estimation for event ", unique(data$Event_ID))
    return(data)
  })
}

#' Perform Parametric Tests
#' @param data Dataframe containing abnormal returns
#' @param event_window Numeric vector specifying event window
#' @param K Number of parameters in market model
#' @return Dataframe with test results
perform_parametric_tests <- function(data, event_window = c(0, 0), K = 4) {
  # Calculate estimation window statistics
  estimation_stats <- data %>%
    filter(Relative_trading_day < 0) %>%
    group_by(Event_ID) %>%
    summarise(
      Mi = sum(!is.na(abnormal_return)),
      S2_ARi = sum(abnormal_return^2, na.rm = TRUE) / (Mi - K),
      S_ARi = sqrt(S2_ARi),
      .groups = 'drop'
    )
  
  # Calculate event window length
  L2 <- event_window[2] - event_window[1] + 1
  
  # Calculate event window returns
  event_stats <- data %>%
    filter(Relative_trading_day >= event_window[1],
           Relative_trading_day <= event_window[2]) %>%
    group_by(Event_ID, Company_name) %>%
    summarise(
      return = if (L2 == 1) {
        mean(abnormal_return, na.rm = TRUE)
      } else {
        sum(abnormal_return, na.rm = TRUE)
      },
      n_obs = sum(!is.na(abnormal_return)),
      .groups = 'drop'
    )
  
  # Calculate test statistics
  results <- event_stats %>%
    left_join(estimation_stats, by = "Event_ID") %>%
    mutate(
      SE = if (L2 == 1) S_ARi else sqrt(L2) * S_ARi,
      t_stat = return / SE,
      df = Mi - K,
      p_value = 2 * (1 - pt(abs(t_stat), df = df)),
      window_type = if (L2 == 1) "AR" else "CAR"
    ) %>%
    mutate(
      significance = case_when(
        p_value < 0.01 ~ "***",
        p_value < 0.05 ~ "**",
        p_value < 0.1  ~ "*",
        TRUE           ~ ""
      )
    )
  
  attr(results, "test_details") <- list(
    window_length = L2,
    test_type = if (L2 == 1) "AR t-test" else "CAR t-test",
    df = "t ~ Mi - K",
    window = paste0("(", event_window[1], ",", event_window[2], ")")
  )
  
  return(results)
}

#' Format Test Results
#' @param results Dataframe with test results
#' @return Formatted dataframe
format_results <- function(results) {
  results %>%
    select(Event_ID, Company_name, return, SE, t_stat, p_value, significance) %>%
    mutate(
      return = round(return * 100, 4),
      SE = round(SE * 100, 4),
      t_stat = round(t_stat, 4),
      p_value = round(p_value, 4)
    ) %>%
    rename(
      `Return (%)` = return,
      `Standard Error (%)` = SE,
      `t-statistic` = t_stat,
      `p-value` = p_value,
      Significance = significance
    )
}

# -----------------------------------------------------------------------------
# Main Execution
# -----------------------------------------------------------------------------
# Create output directories
output_dir <- "output"
output_paths <- list(
  parametric = file.path(output_dir, "parametric")
)
lapply(output_paths, create_dir)

# Read and process data
event_data <- read.csv("dataset")

# Calculate abnormal returns
all_events_data <- event_data %>%
  group_by(Event_ID) %>%
  group_modify(~ calculate_abnormal_returns(.x)) %>%
  ungroup()

# Process results for each window
results_list <- list()
for (window_name in names(EVENT_WINDOWS)) {
  window <- EVENT_WINDOWS[[window_name]]
  
  # Perform tests
  results <- perform_parametric_tests(all_events_data, window, K = 4)
  formatted_results <- format_results(results)
  window_str <- paste0("(", window[1], ",", window[2], ")")
  
  # Save individual CSV
  filename <- paste0("parametric_test_window_", window[1], "_to_", window[2], ".csv")
  write.csv(formatted_results,
            file.path(output_paths$parametric, filename),
            row.names = FALSE)
  
  results_list[[window_str]] <- formatted_results
  
  # Print summary
  cat("\nResults for window", window_str, ":\n")
  cat("Test type:", attr(results, "test_details")$test_type, "\n")
  cat("Window length:", attr(results, "test_details")$window_length, "\n")
  
  sig_summary <- formatted_results %>%
    summarise(
      n_total = n(),
      sig_01 = sum(`p-value` < 0.01, na.rm = TRUE),
      sig_05 = sum(`p-value` < 0.05 & `p-value` >= 0.01, na.rm = TRUE),
      sig_10 = sum(`p-value` < 0.10 & `p-value` >= 0.05, na.rm = TRUE)
    )
  
  cat("Total events:", sig_summary$n_total, "\n")
  cat("Significant at 1%:", sig_summary$sig_01, "\n")
  cat("Significant at 5%:", sig_summary$sig_05, "\n")
  cat("Significant at 10%:", sig_summary$sig_10, "\n\n")
}

# Save combined results
write_xlsx(results_list,
           file.path(output_paths$parametric, "all_parametric_tests.xlsx"))

# Create and save summary table
summary_table <- data.frame(
  Window = names(results_list),
  Events = sapply(results_list, nrow),
  Mean_Return = sapply(results_list, function(x) mean(x$`Return (%)`, na.rm = TRUE)),
  Median_Return = sapply(results_list, function(x) median(x$`Return (%)`, na.rm = TRUE)),
  Sig_1pct = sapply(results_list, function(x) sum(x$`p-value` < 0.01, na.rm = TRUE)),
  Sig_5pct = sapply(results_list, function(x) 
    sum(x$`p-value` < 0.05 & x$`p-value` >= 0.01, na.rm = TRUE)),
  Sig_10pct = sapply(results_list, function(x)
    sum(x$`p-value` < 0.10 & x$`p-value` >= 0.05, na.rm = TRUE))
)

summary_table <- summary_table %>%
  mutate(across(where(is.numeric), ~ round(., 4)))

write.csv(summary_table,
          file.path(output_paths$parametric, "parametric_tests_summary.csv"),
          row.names = FALSE)
