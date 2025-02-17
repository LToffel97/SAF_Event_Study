# Event Study Analysis with Corrado Rank Test
# This script performs event study analysis using the Corrado Rank Test methodology
# to analyze abnormal returns around specific events.

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
  event_data_win  <- data %>% filter(Relative_trading_day >= 0)
  
  # Validate data availability
  if (nrow(estimation_data) == 0 || nrow(event_data_win) == 0) {
    warning("Insufficient data for event ", unique(data$Event_ID))
    return(data)
  }
  
  tryCatch({
    # Fit market model
    market_model <- lm(Stock_return ~ Market_return + World_return + FX_weighted_pct_change,
                      data = estimation_data)
    
    # Calculate predictions and abnormal returns
    estimation_data <- estimation_data %>%
      mutate(
        predicted_return = predict(market_model, newdata = estimation_data),
        abnormal_return  = Stock_return - predicted_return
      )
    
    event_data_win <- event_data_win %>%
      mutate(
        predicted_return = predict(market_model, newdata = event_data_win),
        abnormal_return  = Stock_return - predicted_return
      )
    
    bind_rows(estimation_data, event_data_win)
    
  }, error = function(e) {
    warning("Error in market model estimation for event ", unique(data$Event_ID))
    return(data)
  })
}

#' Perform Corrado Rank Test
#' @param data Dataframe containing abnormal returns
#' @param window Numeric vector specifying event window (start, end)
#' @param output_path Character string specifying output directory
#' @return List containing test results
corrado_rank_test <- function(data, window = c(0, 0), output_path = "output/rank") {
  Mi <- sum(data$Relative_trading_day < 0)  # Estimation window length
  L2 <- window[2] - window[1] + 1           # Event window length
  
  # Calculate scaled ranks
  ranked_data <- data %>%
    group_by(Event_ID) %>%
    mutate(
      rank_AR = rank(abnormal_return),
      K       = rank_AR / (1 + Mi + L2)
    ) %>%
    ungroup()
  
  # Compute daily statistics
  daily_stats <- ranked_data %>%
    filter(Relative_trading_day >= window[1], Relative_trading_day <= window[2]) %>%
    group_by(Relative_trading_day) %>%
    summarise(
      Nt    = n(),
      K_bar = mean(K, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Calculate test statistics
  S_K <- sqrt(sum((daily_stats$K_bar - 0.5)^2) / L2)
  
  daily_stats <- daily_stats %>%
    mutate(
      z_stat       = (K_bar - 0.5) / S_K,
      p_value      = 2 * (1 - pnorm(abs(z_stat))),
      significance = case_when(
        p_value < 0.01 ~ "***",
        p_value < 0.05 ~ "**",
        p_value < 0.1  ~ "*",
        TRUE           ~ ""
      )
    )
  
  # Calculate CAAR statistics
  K_bar_window <- mean(daily_stats$K_bar)
  z_CAAR       <- sqrt(L2) * ((K_bar_window - 0.5) / S_K)
  p_value_CAAR <- 2 * (1 - pnorm(abs(z_CAAR)))
  
  # Prepare results
  results <- list(
    daily_results = daily_stats,
    CAAR_results  = data.frame(
      window       = paste0("(", window[1], ",", window[2], ")"),
      L2           = L2,
      K_bar        = K_bar_window,
      S_K          = S_K,
      z_stat       = z_CAAR,
      p_value      = p_value_CAAR,
      significance = case_when(
        p_value_CAAR < 0.01 ~ "***",
        p_value_CAAR < 0.05 ~ "**",
        p_value_CAAR < 0.1  ~ "*",
        TRUE                 ~ ""
      )
    )
  )
  
  # Save results if output path is provided
  if (!is.null(output_path)) {
    window_str <- paste0("window_", window[1], "_to_", window[2])
    write.csv(daily_stats,
              file = file.path(output_path, paste0("corrado_rank_daily_", window_str, ".csv")),
              row.names = FALSE)
    write.csv(results$CAAR_results,
              file = file.path(output_path, paste0("corrado_rank_CAAR_", window_str, ".csv")),
              row.names = FALSE)
  }
  
  return(results)
}

# -----------------------------------------------------------------------------
# Main Execution
# -----------------------------------------------------------------------------
# Create output directories
output_dir <- "output"
output_paths <- list(
  rank    = file.path(output_dir, "rank"),
  metrics = file.path(output_dir, "metrics")
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
results_list <- lapply(names(EVENT_WINDOWS), function(window_name) {
  cat("\nAnalyzing window", window_name, ":\n")
  corrado_rank_test(
    all_events_data,
    EVENT_WINDOWS[[window_name]],
    output_paths$rank
  )
})
names(results_list) <- names(EVENT_WINDOWS)

# Create and save summary table
summary_table <- do.call(rbind, lapply(names(EVENT_WINDOWS), function(name) {
  results_list[[name]]$CAAR_results %>%
    mutate(window_name = name, .before = 1)
}))

write.csv(
  summary_table,
  file = file.path(output_paths$rank, "corrado_rank_test_summary.csv"),
  row.names = FALSE
)
