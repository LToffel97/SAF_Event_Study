# Cross-sectional T-test for Event Study Analysis
# This script performs cross-sectional t-tests to analyze abnormal returns
# around event dates using a market model approach.

# -----------------------------------------------------------------------------
# Dependencies
# -----------------------------------------------------------------------------
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  dplyr      # Data manipulation
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
    # Estimate market model
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

#' Perform Cross-Sectional Tests
#' @param data Dataframe containing abnormal returns
#' @param window Numeric vector specifying event window (start, end)
#' @param output_path Character string specifying output directory
#' @return List containing test results
perform_cross_sectional_tests <- function(data, window = c(0, 0), output_path = "output/cross_sectional") {
  # Filter event window data
  event_window_data <- data %>%
    filter(Relative_trading_day >= window[1],
           Relative_trading_day <= window[2])
  
  # Compute Daily AAR and Test Statistics
  daily_AARs <- event_window_data %>%
    group_by(Relative_trading_day) %>%
    summarise(
      N      = n(),
      AAR    = mean(abnormal_return, na.rm = TRUE),
      S_AAR  = sd(abnormal_return, na.rm = TRUE),
      t_stat = sqrt(N) * AAR / S_AAR,
      df     = N - 1,
      p_value = 2 * (1 - pt(abs(t_stat), df)),
      .groups = 'drop'
    ) %>%
    mutate(
      significance = case_when(
        p_value < 0.01 ~ "***",
        p_value < 0.05 ~ "**",
        p_value < 0.1  ~ "*",
        TRUE           ~ ""
      )
    )
  
  # Compute CAAR
  event_CARs <- event_window_data %>%
    group_by(Event_ID) %>%
    summarise(CAR = sum(abnormal_return, na.rm = TRUE), .groups = 'drop')

  # Cross-Sectional Test for CAAR
  CAAR_test <- event_CARs %>%
    summarise(
      N      = n(),
      CAAR   = mean(CAR, na.rm = TRUE),
      S_CAAR = sd(CAR, na.rm = TRUE),
      t_stat = sqrt(N) * CAAR / S_CAAR,
      df     = N - 1,
      p_value = 2 * (1 - pt(abs(t_stat), df))
    ) %>%
    mutate(
      significance = case_when(
        p_value < 0.01 ~ "***",
        p_value < 0.05 ~ "**",
        p_value < 0.1  ~ "*",
        TRUE           ~ ""
      )
    )
  
  # Format results
  daily_results <- daily_AARs %>%
    mutate(across(c(AAR, S_AAR), ~ round(. * 100, 4)),
           across(c(t_stat, p_value), ~ round(., 4)))
  
  CAAR_results <- CAAR_test %>%
    mutate(across(c(CAAR, S_CAAR), ~ round(. * 100, 4)),
           across(c(t_stat, p_value), ~ round(., 4)))
  
  # Save results
  if (!is.null(output_path)) {
    window_str <- paste0("window_", window[1], "_to_", window[2])
    write.csv(daily_results,
              file.path(output_path, paste0("cross_sectional_AAR_", window_str, ".csv")),
              row.names = FALSE)
    write.csv(CAAR_results,
              file.path(output_path, paste0("cross_sectional_CAAR_", window_str, ".csv")),
              row.names = FALSE)
  }
  
  return(list(
    daily_results = daily_results,
    CAAR_results  = CAAR_results
  ))
}

# -----------------------------------------------------------------------------
# Main Execution
# -----------------------------------------------------------------------------
# Create output directories
output_dir <- "output"
output_paths <- list(
  cross_sectional = file.path(output_dir, "cross_sectional")
)
lapply(output_paths, create_dir)

# Read and process data
event_data <- read.csv("dataset")

# Calculate abnormal returns
all_events_data <- event_data %>%
  group_by(Event_ID) %>%
  group_modify(~ calculate_abnormal_returns(.x)) %>%
  ungroup()

# Run analysis for each window
results_list <- list()
for (window_name in names(EVENT_WINDOWS)) {
  cat("\nAnalyzing window", window_name, ":\n")
  results_list[[window_name]] <- perform_cross_sectional_tests(
    all_events_data,
    EVENT_WINDOWS[[window_name]],
    output_paths$cross_sectional
  )
}

# Create and save summary table
summary_table <- data.frame(
  Window       = names(EVENT_WINDOWS),
  Window_Size  = sapply(EVENT_WINDOWS, function(w) paste0("(", w[1], ",", w[2], ")")),
  N            = sapply(results_list, function(r) r$CAAR_results$N),
  CAAR         = sapply(results_list, function(r) r$CAAR_results$CAAR),
  t_stat       = sapply(results_list, function(r) r$CAAR_results$t_stat),
  p_value      = sapply(results_list, function(r) r$CAAR_results$p_value),
  significance = sapply(results_list, function(r) r$CAAR_results$significance)
)

write.csv(
  summary_table,
  file.path(output_paths$cross_sectional, "cross_sectional_tests_summary.csv"),
  row.names = FALSE
)
