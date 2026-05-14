# Load required library
library(tidyverse)

# Set the file path 
input_file <- "C:/Users/amrit/Downloads/DEXUSUK (1).csv"
output_file <- "C:/Users/amrit/Downloads/DEXUSUK (1)_cleaned.csv"

# Read the CSV file
data <- read.csv(input_file, stringsAsFactors = FALSE)

# Display number of rows before cleaning
cat("Number of rows before cleaning:", nrow(data), "\n")

# Display number of missing values per column
cat("Missing values per column:\n")
print(colSums(is.na(data)))

# Remove rows with any missing values
clean_data <- na.omit(data)


cat("Number of rows after cleaning:", nrow(clean_data), "\n")


write.csv(clean_data, output_file, row.names = FALSE)

cat("Cleaned data saved to:", output_file, "\n")

# --- Load Libraries ---
library(dplyr)
library(ggplot2)
library(tseries)
library(readr)
library(scales)
library(Metrics)
library(lubridate)
library(zoo)
library(tidyr)
library(forecast)
library(reshape2)
library(rugarch)

# Load Data 
data <- read.csv("C:/Users/amrit/Downloads/DEXUSUK (1)_cleaned.csv",
                 header = TRUE, stringsAsFactors = FALSE)
data$Date <- as.Date(data$Date, format = "%d-%m-%Y")

# Split into training and test datasets
train_data <- filter(data, Date >= as.Date("2015-06-01") & Date < as.Date("2024-06-01"))
test_data  <- filter(data, Date >= as.Date("2024-06-01") & Date <= as.Date("2025-06-30"))

# Log-transform the series
train_data$LogValue <- log(train_data$DEXUSUK)

# Compute first differences
train_data$Diff <- c(NA, diff(train_data$LogValue, differences = 1))

# Estimate slope a 
a <- mean(train_data$Diff, na.rm = TRUE)
cat("Estimated slope a (average drift per period):", a, "\n")

# Plot first differenced series 
ggplot(train_data, aes(x = Date, y = Diff)) +
  geom_line(color = "gray40") +
  geom_hline(yintercept = a, color = "black", linetype = "dashed", size = 0.8) +
  labs(title = "Differenced Log Series",
       y = "Difference") +
  theme_minimal()

# Augmented Dickey-Fuller Test 
adf_result <- adf.test(diff_series, alternative = "stationary")
print(adf_result)

# KPSS Test 
kpss_result <- kpss.test(diff_series, null = "Level")  
print(kpss_result)

# Compute log-returns (epsilon_t) 
epsilon_t <- diff(log(train_data$DEXUSUK))

# Define range 
k_values <- 100:130

# Compute proportion of scaled epsilon_t 
proportions <- numeric(length(k_values))

for (i in seq_along(k_values)) {
  k <- k_values[i]
  epsilon_scaled <- k * epsilon_t
  proportions[i] <- mean(epsilon_scaled >= -1 & epsilon_scaled <= 1)
}

# find optimal K
results_df <- data.frame(K = k_values, Proportion = proportions)
print(results_df)

# Set K 
k <- 125

# Scale and round
epsilon_scaled <- k * epsilon_t
epsilon_discrete <- round(epsilon_scaled)  

# Create Frequency Table and Probabilities 
table_eps <- table(epsilon_discrete)
print(table_eps)

prob_eps <- prop.table(table_eps)
print(prob_eps)


# Original probabilities from your table
prob_eps <- c(
  `-10` = 0.0004476276, `-4` = 0.0008952551, `-3` = 0.0040286482,
  `-2` = 0.0210384960, `-1` = 0.1857654432, `0` = 0.5778871979,
  `1` = 0.1866606983, `2` = 0.0192479857, `3` = 0.0035810206,
  `4` = 0.0004476276
)

# Define support [-2, 2]
core_vals <- as.character(-2:2)
core_mass <- prob_eps[core_vals]
core_mass[is.na(core_mass)] <- 0
total_core_mass <- sum(core_mass)

# Calculate excess mass from outside [-2, 2]
outside_vals <- setdiff(names(prob_eps), core_vals)
outside_mass <- sum(prob_eps[outside_vals])

# Redistribute excess mass proportionally inside [-2, 2]
adjusted_core_mass <- core_mass + (core_mass / total_core_mass) * outside_mass
adjusted_core_mass <- adjusted_core_mass / sum(adjusted_core_mass)  # Normalize to 1

# Create data frame for plotting
df_plot <- data.frame(
  Epsilon = factor(core_vals, levels = core_vals),
  Probability = as.numeric(adjusted_core_mass)
)

# Original proportions from epsilon_t
round(prob_eps[support_vals_str] * 100, 2)

# Plot histogram 
ggplot(df_plot, aes(x = Epsilon, y = Probability)) +
  geom_bar(stat = "identity", fill = "steelblue", color = "black") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(title = "εt Distribution",
       x = "εt Value", y = "Percentage") +
  theme_minimal()

# Adjusted probabilities from your earlier step
adjusted_core_mass <- c(
  `-2` = 0.021238,
  `-1` = 0.187528,
  `0`  = 0.583371,
  `1`  = 0.188432,
  `2`  = 0.019431
)

# Calculate q and p as averages of symmetric probabilities
q <- (adjusted_core_mass["-2"] + adjusted_core_mass["2"]) / 2
p <- (adjusted_core_mass["-1"] + adjusted_core_mass["1"]) / 2

# Calculate the remaining probability mass for 0
p0 <- 1 - 2 * p - 2 * q

# Print results
cat("q =", round(q, 6), "\n")
cat("p =", round(p, 6), "\n")
cat("1 - 2p - 2q =", round(p0, 6), "\n")

# Define the support and corresponding probabilities
epsilon_values <- c(-2, -1, 0, 1, 2)

# Use the symmetric probabilities based on earlier p and q
q <- 0.020334
p <- 0.18798
p0 <- 1 - 2 * p - 2 * q  

epsilon_probs <- c(q, p, p0, p, q)

# Generate  samples
epsilon_sim <- sample(epsilon_values, size = 1000, replace = TRUE, prob = epsilon_probs)

# Set up discrete distribution
epsilon_values <- c(-2, -1, 0, 1, 2)
q <- 0.020334
p <- 0.18798
p0 <- 1 - 2 * p - 2 * q  

epsilon_probs <- c(q, p, p0, p, q)
names(epsilon_probs) <- as.character(epsilon_values)  

round(epsilon_probs[support_vals_str] * 100, 2)


# Simulation function
simulate_exchange_rate_log <- function(T, N_paths = 1000, S0 = 1.0) {
  eps_matrix <- matrix(sample(epsilon_values, size = T * N_paths,
                              replace = TRUE, prob = epsilon_probs),
                       nrow = T, ncol = N_paths)
  # Scale epsilon
  eps_matrix <- eps_matrix / k  
  
  log_paths <- apply(eps_matrix, 2, cumsum)
  log_paths <- log_paths + log(S0)
  return(log_paths)
}

# Short-term plotting function
plot_short_log_paths <- function(log_matrix, title_str, actual_log = NULL) {
  T <- nrow(log_matrix)
  N <- ncol(log_matrix)
  days <- 1:T
  
  sample_indices <- sample(N, 5)
  sampled_paths <- log_matrix[, sample_indices]
  
  path_df <- as.data.frame(sampled_paths)
  path_df$Day <- days
  path_df_long <- pivot_longer(path_df, -Day, names_to = "Path", values_to = "Log_Rate")
  
  plt <- ggplot() +
    geom_line(data = path_df_long, aes(x = Day, y = Log_Rate, group = Path),
              color = "darkred", size = 0.8, alpha = 0.7) +
    scale_x_continuous(breaks = 1:T) +  # integer-only x-axis for short-term
    labs(title = title_str, x = "Day", y = "Log Exchange Rate") +
    theme_minimal(base_size = 13) +
    theme(legend.position = "none",
          plot.title = element_text(size = 15, face = "bold"))
  
  if(!is.null(actual_log)) {
    df_actual <- data.frame(Day = 1:length(actual_log), Log_Rate = actual_log)
    plt <- plt + geom_line(data = df_actual, aes(x = Day, y = Log_Rate),
                           color = "black", size = 1.0)
  }
  
  return(plt)
}

# Plotting function 
plot_log_paths <- function(log_matrix, title_str, actual_log = NULL) {
  T <- nrow(log_matrix)
  N <- ncol(log_matrix)
  days <- 1:T
  
  sample_indices <- sample(N, 5)
  sampled_paths <- log_matrix[, sample_indices]
  
  path_df <- as.data.frame(sampled_paths)
  path_df$Day <- days
  path_df_long <- pivot_longer(path_df, -Day, names_to = "Path", values_to = "Log_Rate")
  
  plt <- ggplot() +
    geom_line(data = path_df_long, aes(x = Day, y = Log_Rate, group = Path),
              color = "darkred", size = 0.8, alpha = 0.7) +
    labs(title = title_str, x = "Day", y = "Log Exchange Rate") +
    theme_minimal(base_size = 13) +
    theme(legend.position = "none",
          plot.title = element_text(size = 15, face = "bold"))
  
  if(!is.null(actual_log)) {
    df_actual <- data.frame(Day = 1:length(actual_log), Log_Rate = actual_log)
    plt <- plt + geom_line(data = df_actual, aes(x = Day, y = Log_Rate),
                           color = "black", size = 1.0)
  }
  
  return(plt)
}

# Define test data 
log_test_short <- log(test_data %>% slice(1:3)   %>% pull(DEXUSUK))
log_test_mid   <- log(test_data %>% slice(1:90)  %>% pull(DEXUSUK))
log_test_long  <- log(test_data %>% slice(1:250) %>% pull(DEXUSUK))

S0 <- test_data$DEXUSUK[1]   

# Simulate log paths
short_log_sim <- simulate_exchange_rate_log(T = 3,   N_paths = 1000, S0 = S0)
mid_log_sim   <- simulate_exchange_rate_log(T = 90,  N_paths = 1000, S0 = S0)
long_log_sim  <- simulate_exchange_rate_log(T = 250, N_paths = 1000, S0 = S0)

# Step 8: Plot simulated paths with actual test data
plot_short <- plot_short_log_paths(short_log_sim, "Short-Term (3 Days)", actual_log = log_test_short)
plot_mid   <- plot_log_paths(mid_log_sim, "Mid-Term (90 Days)", actual_log = log_test_mid)
plot_long  <- plot_log_paths(long_log_sim, "Long-Term (250 Days)", actual_log = log_test_long)

print(plot_short)
print(plot_mid)
print(plot_long)

# Compute RMSE in log scale ---
short_log_mean <- rowMeans(short_log_sim)
mid_log_mean   <- rowMeans(mid_log_sim)
long_log_mean  <- rowMeans(long_log_sim)

rmse_short_log <- rmse(log_test_short, short_log_mean)
rmse_mid_log   <- rmse(log_test_mid,   mid_log_mean)
rmse_long_log  <- rmse(log_test_long,  long_log_mean)

cat("RMSE Short-term (log, DEXUSUK):", rmse_short_log, "\n")
cat("RMSE Mid-term (log, DEXUSUK):",   rmse_mid_log, "\n")
cat("RMSE Long-term (log, DEXUSUK):",  rmse_long_log, "\n")

# RMSE in % terms
cat("RMSE Short-term (log %):", rmse_short_log * 100, "%\n")
cat("RMSE Mid-term (log %):",   rmse_mid_log * 100, "%\n")
cat("RMSE Long-term (log %):",  rmse_long_log * 100, "%\n")

log_x <- train_data$LogValue
diff_logx <- na.omit(diff(log_x))

# Auto ARIMA 
auto_fit <- auto.arima(log_x, seasonal = FALSE, stepwise = FALSE, approximation = FALSE)
summary(auto_fit)

par(mfrow=c(1,2))
acf(diff_logx, main="ACF of Differenced Series")   
pacf(diff_logx, main="PACF of Differenced Series") 


# Simulations 
simulate_arima_rw_log <- function(fit, T, N_paths = 1000, S0_log = 0) {
  sim_matrix <- matrix(NA, nrow = T, ncol = N_paths)
  for (i in 1:N_paths) {
    sim_log <- simulate(fit, nsim = T, future = TRUE)
    # Align starting value with log(S0)
    sim_log <- sim_log - sim_log[1] + S0_log
    sim_matrix[, i] <- sim_log
  }
  return(sim_matrix)
}

plot_fan_chart_log <- function(sim_result, title, actual_log) {
  days <- seq_len(nrow(sim_result))
  
  
  idx <- sample(ncol(sim_result), 5)
  sim_subset <- sim_result[, idx, drop = FALSE]
  df_paths <- reshape2::melt(sim_subset)
  colnames(df_paths) <- c("Day", "Path", "Log_Value")
  
  # Actual log test data
  df_actual <- data.frame(Day = seq_along(actual_log), Log_Value = actual_log)
  
  p <- ggplot() +
    geom_line(data = df_paths, aes(x = Day, y = Log_Value, group = Path),
              color = "red", alpha = 0.6, size = 0.8) +
    geom_line(data = df_actual, aes(x = Day, y = Log_Value),
              color = "black", size = 1.0) +
    labs(title = title, x = "Day", y = "Log Exchange Rate (DEXUSUK)") +
    theme_minimal()
  
  return(p)
}

# Horizons
horizons <- c(3, 90, 250)

# Initialize lists
simulations <- list()
rmse_vals <- numeric(length(horizons))
rmse_pct  <- numeric(length(horizons))   
plots     <- list()

# Run simulation, compute RMSE, and plot
for (i in seq_along(horizons)) {
  T_h <- horizons[i]
  
  # Simulate log paths (S0_log = log of initial DEXUSUK)
  sim_log <- simulate_arima_rw_log(fit_log, T = T_h, N_paths = 1000, S0_log = log(test_data$DEXUSUK[1]))
  simulations[[i]] <- sim_log
  
  # Select actual log test data
  actual_log <- if(T_h == 3) log_test_short else if(T_h == 90) log_test_mid else log_test_long
  
  # RMSE in log scale
  rmse_vals[i] <- rmse(actual_log, rowMeans(sim_log))
  
  # Rescaled RMSE in %
  rmse_pct[i] <- rmse_vals[i] * 100
  
  # Create fan chart plot
  p <- plot_fan_chart_log(sim_log,
                          title = paste0("Auto ARIMA Log Simulation - ", T_h, " Days"),
                          actual_log = actual_log)
  
  # For short-term (3 days), fix integer x-axis
  if(T_h == 3) {
    p <- p + scale_x_continuous(breaks = 1:3)
  }
  
  plots[[i]] <- p
}

# Name RMSE values
names(rmse_vals) <- names(rmse_pct) <- c("Short-term (3 days)", "Mid-term (90 days)", "Long-term (250 days)")

# Print RMSE in log units and rescaled %
cat("RMSE (log scale):\n")
print(rmse_vals)
cat("\nRMSE (log scale, %):\n")
print(rmse_pct)

# Show plots
print(plots[[1]])
print(plots[[2]])
print(plots[[3]])


# Log returns
log_returns <- diff(log(train_data$DEXUSUK))
mu_hat <- mean(log_returns)
sigma_hat <- sd(log_returns)

cat("Estimated mu:", mu_hat, "\n")
cat("Estimated sigma:", sigma_hat, "\n")

# Simulation function 
simulate_white_noise_rw_log <- function(T, N_paths = 1000, S0 = 1.0, mu = -9.70792e-05, sigma = 0.006094623) {
  eps_matrix <- matrix(rnorm(T * N_paths, mean = mu, sd = sigma),
                       nrow = T, ncol = N_paths)
  
  log_paths <- apply(eps_matrix, 2, cumsum)
  log_paths <- log_paths + log(S0)
  
  return(log_paths)  
}

# plotting function 
plot_fan_chart_log <- function(log_matrix, title_str, actual_log, short_term = FALSE) {
  T <- nrow(log_matrix)
  N <- ncol(log_matrix)
  days <- 1:T
  
  # Sample 5 paths for display
  sample_indices <- sample(N, 5)
  sampled_paths <- log_matrix[, sample_indices]
  
  path_df <- as.data.frame(sampled_paths)
  path_df$Day <- days
  path_df_long <- pivot_longer(path_df, -Day, names_to = "Path", values_to = "Log_Rate")
  
  # Actual test data
  actual_df <- data.frame(Day = 1:length(actual_log), Log_Rate = actual_log)
  
  plt <- ggplot(path_df_long, aes(x = Day, y = Log_Rate, group = Path)) +
    geom_line(size = 0.8, color = "steelblue", alpha = 0.7) +
    geom_line(data = actual_df, aes(x = Day, y = Log_Rate),
              inherit.aes = FALSE, color = "black", size = 1.0) +
    labs(title = title_str, x = "Day", y = "Log Exchange Rate") +
    theme_minimal(base_size = 13) +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 15, face = "bold")
    )
  
  if(short_term) plt <- plt + scale_x_continuous(breaks = 1:T)
  
  return(plt)
}

# Simulate for short, mid, long horizons 
S0 <- test_data$DEXUSUK[1]

short_rw_log <- simulate_white_noise_rw_log(T = 3, N_paths = 1000, S0 = S0, mu = mu_hat, sigma = sigma_hat)
mid_rw_log   <- simulate_white_noise_rw_log(T = 90, N_paths = 1000, S0 = S0, mu = mu_hat, sigma = sigma_hat)
long_rw_log  <- simulate_white_noise_rw_log(T = 250, N_paths = 1000, S0 = S0, mu = mu_hat, sigma = sigma_hat)

# Define test data
test_short_log <- log(test_data$DEXUSUK[1:3])
test_mid_log   <- log(test_data$DEXUSUK[1:90])
test_long_log  <- log(test_data$DEXUSUK[1:250])

# Plot charts 
plot_short_rw <- plot_fan_chart_log(short_rw_log, "White Noise RW - Short-Term (3 Days)", test_short_log, short_term = TRUE)
plot_mid_rw   <- plot_fan_chart_log(mid_rw_log, "White Noise RW - Mid-Term (90 Days)", test_mid_log)
plot_long_rw  <- plot_fan_chart_log(long_rw_log, "White Noise RW - Long-Term (250 Days)", test_long_log)

print(plot_short_rw)
print(plot_mid_rw)
print(plot_long_rw)

# Compute RMSE 
rmse_short_rw <- rmse(test_short_log, rowMeans(short_rw_log))
rmse_mid_rw   <- rmse(test_mid_log, rowMeans(mid_rw_log))
rmse_long_rw  <- rmse(test_long_log, rowMeans(long_rw_log))

cat("White Noise RW RMSE (log) - Short-term:", rmse_short_rw, "\n")
cat("White Noise RW RMSE (log) - Mid-term:", rmse_mid_rw, "\n")
cat("White Noise RW RMSE (log) - Long-term:", rmse_long_rw, "\n")

# RMSE in %
cat("White Noise RW RMSE (log %) - Short-term:", rmse_short_rw*100, "%\n")
cat("White Noise RW RMSE (log %) - Mid-term:", rmse_mid_rw*100, "%\n")
cat("White Noise RW RMSE (log %) - Long-term:", rmse_long_rw*100, "%\n")


# Step 1: Compute log-returns
log_returns <- diff(log(train_data$DEXUSUK))
log_returns <- na.omit(log_returns)

# Information criteria 
infocriteria(fit11)
infocriteria(fit12)

# Specify GARCH(1,1) model 
garch_spec <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
  mean.model = list(armaOrder = c(0, 0), include.mean = TRUE),
  distribution.model = "norm"
)

# Fit GARCH model
garch_fit <- ugarchfit(spec = garch_spec, data = log_returns)
print(garch_fit)

# Define simulation function 
simulate_garch_paths <- function(T, N_paths = 1000, S0 = 1.0, garch_fit) {
  sim_matrix <- matrix(NA, nrow = T, ncol = N_paths)
  
  for (i in 1:N_paths) {
    sim <- ugarchsim(garch_fit, n.sim = T, m.sim = 1)
    sim_log_returns <- fitted(sim)
    sim_matrix[, i] <- cumsum(sim_log_returns) + log(S0)
  }
  
  return(sim_matrix)
}

# Fan chart plotting function
plot_fan_chart_log <- function(log_matrix, title_str, actual_log = NULL, short_term = FALSE) {
  T <- nrow(log_matrix)
  N <- ncol(log_matrix)
  days <- 1:T
  
  sample_indices <- sample(N, 5)
  sampled_paths <- log_matrix[, sample_indices]
  
  path_df <- as.data.frame(sampled_paths)
  path_df$Day <- days
  path_df_long <- pivot_longer(path_df, -Day, names_to = "Path", values_to = "Log_Rate")
  
  plt <- ggplot(path_df_long, aes(x = Day, y = Log_Rate, group = Path)) +
    geom_line(size = 0.8, color = "orange", alpha = 0.7)
  
  if(!is.null(actual_log)) {
    df_actual <- data.frame(Day = 1:length(actual_log), Log_Rate = actual_log)
    plt <- plt + geom_line(data = df_actual, aes(x = Day, y = Log_Rate),
                           inherit.aes = FALSE, color = "black", size = 1.0)
  }
  
  plt <- plt + labs(title = title_str, x = "Day", y = "Log Exchange Rate") +
    theme_minimal(base_size = 13) +
    theme(legend.position = "none",
          plot.title = element_text(size = 15, face = "bold"))
  
  if(short_term) plt <- plt + scale_x_continuous(breaks = 1:T)
  
  return(plt)
}

# Define horizons and test sets 
horizons <- c(3, 90, 250)
S0 <- test_data$DEXUSUK[1]

log_test_short <- log(test_data$DEXUSUK[1:3])
log_test_mid   <- log(test_data$DEXUSUK[1:90])
log_test_long  <- log(test_data$DEXUSUK[1:250])

test_sets <- list(log_test_short, log_test_mid, log_test_long)
names(test_sets) <- c("Short-term (3 days)", "Mid-term (90 days)", "Long-term (250 days)")

# Simulate GARCH paths and compute RMSE ---
garch_simulations <- list()
rmse_garch <- c()
rmse_garch_pct <- c()
plots_garch <- list()

for(i in seq_along(horizons)) {
  T_h <- horizons[i]
  
  sim_log <- simulate_garch_paths(T = T_h, N_paths = 1000, S0 = S0, garch_fit = garch_fit)
  garch_simulations[[i]] <- sim_log
  
  actual <- test_sets[[i]]
  
  # RMSE in log scale
  rmse_garch[i] <- rmse(actual, rowMeans(sim_log))
  # Rescaled to %
  rmse_garch_pct[i] <- rmse_garch[i] * 100
  
  # Fan chart plot
  plots_garch[[i]] <- plot_fan_chart_log(sim_log,
                                         title_str = paste0("GARCH(1,1) Simulation - ", T_h, " Days"),
                                         actual_log = actual,
                                         short_term = (T_h == 3))
}

names(rmse_garch) <- names(rmse_garch_pct) <- names(test_sets)

# Print RMSE 
cat("GARCH(1,1) RMSE (log scale):\n")
print(rmse_garch)

cat("\nGARCH(1,1) RMSE (rescaled %):\n")
print(rmse_garch_pct)

# Show plots
print(plots_garch[[1]])  
print(plots_garch[[2]])  
print(plots_garch[[3]])  

# Standardized residuals from GARCH fit
std_resid <- residuals(garch_fit, standardize = TRUE)

# Ljung-Box test
Box.test(std_resid, lag = 20, type = "Ljung-Box")

# ARCH LM test on standardized residuals
library(FinTS)
ArchTest(std_resid, lags = 20)

