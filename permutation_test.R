# Permutation Test for Event Study Analysis
# This script performs permutation tests to analyze abnormal returns
# and cumulative abnormal returns (CAR) around event dates.

# -----------------------------------------------------------------------------
# Dependencies
# -----------------------------------------------------------------------------
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  dplyr,      # Data manipulation
  tidyr       # Data tidying
)

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
# Set default options
options(stringsAsFactors = FALSE)
options(scipen = 999)

# Constants
N_PERMUTATIONS <- 10000
K_PARAMETERS <- 4  # Number of parameters in market model

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

#' Calculate Market Model for Event Window
#' @param data Dataframe containing event data
#' @param est_start Start of estimation window
#' @param est_end End of estimation window
#' @param event_start Start of event window
#' @param event_end End of event window
#' @return List containing model results
calculate_market_model <- function(data, est_start, est_end, event_start, event_end) {
  # Filter estimation window data
  est_data <- data %>% 
    filter(Relative_trading_day <= est_end,
           Relative_trading_day >= est_start) %>%
    drop_na()
  
  # Filter event window data
  event_data <- data %>% 
    filter(Relative_trading_day <= event_end,
           Relative_trading_day >= event_start)
  
  # Fit market model
  market_model <- lm(Stock_return ~ Market_return + World_return + FX_weighted_pct_change,
                     data = est_data)
  
  # Calculate predictions and abnormal returns
  predicted_returns_est <- predict(market_model, newdata = est_data)
  predicted_returns_event <- predict(market_model, newdata = event_data)
  
  ar_est <- est_data$Stock_return - predicted_returns_est
  ar_event <- event_data$Stock_return - predicted_returns_event
  
  list(
    ar_est = ar_est,
    ar_event = ar_event,
    est_data = est_data,
    event_data = event_data
  )
}

#' Calculate Test Statistics
#' @param ar_data Abnormal returns data
#' @param window_size Size of event window
#' @param k Number of parameters in market model
#' @return Named vector of test statistics
calculate_test_statistics <- function(ar_data, window_size, k = K_PARAMETERS) {
  m_val <- length(ar_data)
  df <- m_val - k
  
  # Calculate variance
  var_ar_est <- (1 / df) * sum(ar_data^2)
  std_ar_est <- sqrt(var_ar_est)
  
  # Adjust variance for window size
  var_car_est <- window_size * var_ar_est
  std_car_est <- sqrt(var_car_est)
  
  c(var = var_ar_est, std = std_ar_est, var_car = var_car_est, std_car = std_car_est)
}

#' Perform Permutation Test
#' @param data Combined abnormal returns data
#' @param observed_stat Observed test statistic
#' @param window_size Size of event window
#' @param n_perm Number of permutations
#' @return List containing test results
perform_permutation_test <- function(data, observed_stat, window_size, n_perm = N_PERMUTATIONS) {
  perm_stats <- numeric(n_perm)
  
  for (i in 1:n_perm) {
    # Shuffle data
    shuffled_data <- sample(data, replace = FALSE)
    
    # Calculate test statistic for permuted data
    stats <- calculate_test_statistics(shuffled_data[-1], window_size)
    ar_perm <- shuffled_data[1]
    
    if (window_size == 1) {
      perm_stats[i] <- ar_perm / stats["std"]
    } else {
      perm_stats[i] <- sum(shuffled_data[1:window_size]) / stats["std_car"]
    }
  }
  
  # Calculate p-value
  p_value <- sum(abs(perm_stats) >= abs(observed_stat), na.rm = TRUE) / n_perm
  
  list(
    perm_stats = perm_stats,
    p_value = p_value,
    observed_stat = observed_stat
  )
}

# -----------------------------------------------------------------------------
# Main Execution
# -----------------------------------------------------------------------------
# Create output directory
output_dir <- "output"
output_path <- file.path(output_dir, "permutation")
create_dir(output_path)

# Read data
df <- read.csv("dataset")

# Initialize results dataframe
results_df <- data.frame(
  Event_ID = character(),
  Window = character(),
  T_obs = numeric(),
  T_perm_median = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

# Process each event
event_ids <- unique(df$Event_ID)

for (id in event_ids) {
  cat("\nProcessing event", id, "...\n")
  
  # Define windows
  windows <- list(
    AR = list(
      est_start = -121, est_end = -1,
      event_start = 0, event_end = 0,
      size = 1
    ),
    CAR_1 = list(
      est_start = -122, est_end = -2,
      event_start = -1, event_end = 1,
      size = 3
    ),
    CAR_2 = list(
      est_start = -123, est_end = -3,
      event_start = -2, event_end = 2,
      size = 5
    )
  )
  
  # Process each window type
  for (win_name in names(windows)) {
    win <- windows[[win_name]]
    
    # Calculate market model
    model_results <- calculate_market_model(
      df[df$Event_ID == id, ],
      win$est_start, win$est_end,
      win$event_start, win$event_end
    )
    
    # Calculate observed statistics
    stats <- calculate_test_statistics(model_results$ar_est, win$size)
    
    if (win$size == 1) {
      observed_stat <- model_results$ar_event[1] / stats["std"]
    } else {
      observed_stat <- sum(model_results$ar_event) / stats["std_car"]
    }
    
    # Perform permutation test
    all_ar <- c(model_results$ar_est, model_results$ar_event)
    perm_results <- perform_permutation_test(all_ar, observed_stat, win$size)
    
    # Store results
    results_df <- rbind(results_df, data.frame(
      Event_ID = id,
      Window = win_name,
      T_obs = perm_results$observed_stat,
      T_perm_median = median(perm_results$perm_stats, na.rm = TRUE),
      p_value = perm_results$p_value,
      stringsAsFactors = FALSE
    ))
  }
}

# Save results by window type
results_by_window <- split(results_df, results_df$Window)
for (win in names(results_by_window)) {
  filename <- paste0("permutation_test_", win, ".csv")
  write.csv(
    results_by_window[[win]],
    file.path(output_path, filename),
    row.names = FALSE
  )
}

# Create and save summary table
summary_table <- results_df %>%
  group_by(Window) %>%
  summarise(
    N = n(),
    Mean_T_stat = mean(T_obs, na.rm = TRUE),
    Median_T_stat = median(T_obs, na.rm = TRUE),
    Sig_1pct = sum(p_value < 0.01, na.rm = TRUE),
    Sig_5pct = sum(p_value < 0.05 & p_value >= 0.01, na.rm = TRUE),
    Sig_10pct = sum(p_value < 0.10 & p_value >= 0.05, na.rm = TRUE)
  ) %>%
  mutate(across(where(is.numeric), ~ round(., 4)))

write.csv(
  summary_table,
  file.path(output_path, "permutation_test_summary.csv"),
  row.names = FALSE
)
