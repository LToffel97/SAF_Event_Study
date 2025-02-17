# SAF_Event_Study

This repository contains a collection of R scripts for conducting event study analysis using different statistical approaches. The implementation includes various test statistics commonly used in event study methodology, including Corrado Rank Test, Cross-sectional Test, Parametric T-test, and Permutation Test.

### Scripts

1. Corrado Rank Test (corrado_rank_test.R)
Implements the non-parametric Corrado Rank Test methodology:

- Calculates abnormal returns using a four-factor market model
- Performs rank testing across multiple event windows
- Generates daily statistics and cumulative average abnormal returns (CAAR)

2. Cross-sectional Test (cross_sectional_test.R)
Performs cross-sectional analysis of abnormal returns:
- Computes daily Average Abnormal Returns (AAR) with t-statistics
- Analyzes Cumulative Average Abnormal Returns (CAAR)
- Tests significance across multiple event windows

3. Parametric T-test (parametric_t_test.R)
Implements standard parametric tests for event study analysis:
- Calculates individual event statistics with proper degrees of freedom
- Performs both AR (single day) and CAR (multiple day) tests
- Accounts for estimation window variance in test statistics

4. Permutation Test (permutation_test.R)
Conducts permutation-based statistical tests:
- Implements permutation testing with 10,000 random permutations
- Calculates test statistics for single-day AR and multi-day CAR
- Computes empirical p-values based on permutation distributions

### Dataset
The analysis requires a dataset (dataset.csv) with the following structure:
- Event identifiers and company information
- Market and stock return data
- Event window specifications
- Exchange rate information

Key variables include:
- Event_ID, Company_name, Ticker
- Stock_return, Market_return, World_return
- Relative_trading_day
- FX_weighted_pct_change

### Usage
Place your event study data in a file named dataset.csv in the root directory
Run the desired test script(s)
Results will be saved in the output directory, organized by test type

Each script creates its own subdirectory in the output folder:

output/
├── rank/             # Corrado Rank Test results
├── cross_sectional/  # Cross-sectional Test results
├── parametric/       # Parametric T-test results
└── permutation/      # Permutation Test results

### Requirements

Required R packages:
install.packages(c("dplyr", "tidyr", "writexl"))

### Output Format

Each test generates:
- Individual test results for each event window
- Summary statistics with significance levels
- Combined results in CSV/Excel format

Results include:
- Test statistics
- p-values
- Significance levels (*, **, ***)
- Sample sizes and window specifications

### Acknowledgments

Code cleaning and documentation was supported by Anthropic's Claude AI assistant
Original analysis scripts were developed for academic research
Implementation based on established event study methodology literature

