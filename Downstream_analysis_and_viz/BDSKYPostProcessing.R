library(s20x)
library(boa)
library(Hmisc)
library(miscTools)
library(viridis)
library(dplyr)

# Function to read configuration file
read_config <- function(config_file) {
  config <- readLines(config_file)
  config_list <- list()
  
  log_files <- c()
  log_files_start <- FALSE
  
  for (line in config) {
    # Skip empty lines
    if (trimws(line) == "") next
    
    if (log_files_start) {
      # Collect log file paths until the end of the file or until another key is encountered
      if (grepl(":", line)) {
        # Detected a new key, stop collecting log files
        log_files_start <- FALSE
        # Process the new key-value pair
        key_value <- strsplit(line, ":", fixed=TRUE)[[1]]
        key <- trimws(key_value[1])
        value <- trimws(paste(key_value[-1], collapse=":"))
        config_list[[key]] <- value
      } else {
        # Collect log file path
        log_files <- c(log_files, trimws(line))
        next  # Continue collecting log files
      }
    }
    
    if (!log_files_start) {
      # Check if line indicates start of log_files
      if (trimws(line) == "log_files:" || trimws(line) == "log_files") {
        log_files_start <- TRUE
      } else {
        # Process key-value pairs
        key_value <- strsplit(line, ":", fixed=TRUE)[[1]]
        if (length(key_value) >= 2) {
          key <- trimws(key_value[1])
          value <- trimws(paste(key_value[-1], collapse=":"))
          config_list[[key]] <- value
        } else {
          cat("Invalid line in config file (missing ':'):\n", line, "\n")
        }
      }
    }
  }
  
  if (length(log_files) == 0) {
    stop("No log files specified in the configuration file.")
  }
  
  config_list$log_files <- log_files
  return(config_list)
}

# Get the configuration file path from the user
cat("\nPlease enter the location of the configuration file: ")
config_file <- readLines(file("stdin"), 1)

# Read the configuration file
config <- read_config(config_file)

# Extract parameters from the configuration
burnin_percent <- as.integer(config$burnin_percent)
recent <- as.numeric(config$recent_sample_date)
timechange_array <- as.numeric(unlist(strsplit(config$time_steps, " ")))
log_files <- config$log_files

closeAllConnections()

# Prepare to plot all lineages on one plot
plot(1, ylab=expression(R[e] ~ "(Effective Reproduction Number)"), xlim=range(c(recent - timechange_array, recent)), ylim=c(0,5), 
     xlab="Year", col='white', cex.lab=1.2, cex.axis=1.1)
minor.tick(nx=2, ny=2, tick.ratio=.2)
abline(h=1, col='grey')

# Initialize lists to store median Re values and HPD intervals at each time point for each population
medians_by_population <- list()
hpd_by_population <- list()
F_times <- NULL  # To store the time points (same for all populations)

# Dynamically determine the names based on the number of log files processed
pop_names <- c("NY10", "WN02")  # Adjust or extend as necessary
if (length(log_files) > length(pop_names)) {
  pop_names <- c(pop_names, paste("Pop", (length(pop_names)+1):length(log_files), sep=""))
}

# Define the colors for the plot
first_series_color <- "#7C1158"  # Specific color for the first series
second_series_color <- "#11A579"  # Contrasting color for the second series
colors <- c(first_series_color, second_series_color)  

# Function to process a single log file and plot its data with shaded uncertainty region
process_and_plot_log_file <- function(log_file, i, burnin_percent, timechange_array, recent, colors) {
  # Read the log file
  log_data <- read.table(log_file, header=T)
  
  # Extract column names for reproductive number and related data
  R0_names <- names(log_data)[which(regexpr("reproductiveNumber.", names(log_data)) > 0)]
  
  nsamples <- nrow(log_data)
  burnin <- round(burnin_percent * nsamples / 100)
  
  intervalNumber <- length(timechange_array) + 1  # +1 to account for the recent value
  
  F <- matrix(data=NA, nrow=nsamples - burnin, ncol=intervalNumber)
  
  # Time steps for plotting (forward in time)
  F_times <- c(recent, recent - timechange_array)  # Time axis moves forward from the recent date
  
  # Populate the data matrix (with correct burn-in) for reproductive numbers
  for(k in 1:(nsamples - burnin)) {
    for (l in 1:(intervalNumber-1)) {
      if (l <= length(R0_names)) {
        F[k, l] <- log_data[[R0_names[l]]][k + burnin]
      } else {
        F[k, l] <- log_data[[R0_names[length(R0_names)]]][k + burnin]
      }
    }
    F[k, intervalNumber] <- log_data[[R0_names[length(R0_names)]]][k + burnin]
  }
  
  # Reverse both the time axis and the data
  F <- F[, ncol(F):1]  # Flip the data to match the reversed time
  
  # Calculate medians and HPD intervals for plotting
  medians <- rep(NA, intervalNumber)
  hpd_F <- matrix(data=NA, nrow=2, ncol=intervalNumber)
  
  for(j in 1:intervalNumber) {
    medians[j] <- median(F[, j], na.rm=T)
    hpd_F[, j] <- boa.hpd(F[, j], 0.05)[1:2]  # 95% HPD interval
  }
  
  # Plot the shaded region for the HPD intervals (uncertainty)
  polygon(c(F_times, rev(F_times)), c(hpd_F[2, ], rev(hpd_F[1, ])), col=adjustcolor(colors[i], alpha.f=0.2), border=NA)
  
  # Plot the median line and points
  lines(F_times, medians, col=colors[i], lwd=2)
  points(F_times, medians, pch=16, cex=1.2, col=colors[i])
  
  # Return the median Re values, HPD intervals, and time points
  return(list(medians = medians, hpd_F = hpd_F, F_times = F_times))
}

# Loop through log files and process each one
hpd_by_population <- list()
for(i in 1:length(log_files)) {
  res <- process_and_plot_log_file(log_files[i], i, burnin_percent, timechange_array, recent, colors)
  
  medians_by_population[[i]] <- res$medians  # Store the median Re values at each time point for this population
  hpd_by_population[[i]] <- res$hpd_F       # Store the HPD intervals
  
  if (is.null(F_times)) {
    F_times <- res$F_times  # Store the time points (same for all populations)
  }
}

# Add a legend to identify lineages
legend("topright", legend=pop_names[1:length(log_files)], col=colors[1:length(log_files)], lwd=2, pch=16, bty="n")

# Separate summer and winter indices
num_time_points <- length(F_times)
summer_indices <- seq(1, num_time_points, by=2)
winter_indices <- seq(2, num_time_points, by=2)

# Initialize lists to store medians for boxplots
median_Re_summer <- list()
median_Re_winter <- list()

for (i in 1:length(log_files)) {
  medians <- medians_by_population[[i]]
  
  median_Re_summer[[i]] <- medians[summer_indices]
  median_Re_winter[[i]] <- medians[winter_indices]
}

# Prepare data for boxplots
boxplot_data <- list(
  "NY10 Summer" = median_Re_summer[[1]],
  "NY10 Winter" = median_Re_winter[[1]],
  "WN02 Summer" = median_Re_summer[[2]],
  "WN02 Winter" = median_Re_winter[[2]]
)

# Assign colors accordingly
boxplot_colors <- c(colors[1], colors[1], colors[2], colors[2])

# Create the boxplot
par(mfrow=c(1,1))
boxplot(boxplot_data,
        main="Median Re Values Across Time Points",
        col=boxplot_colors,
        ylab=expression(R[e]),
        names=names(boxplot_data),
        las=2)

# --- calculate and print some statistics ---

# Initialize vectors to store all data
Re_values <- c()
Population <- c()
Season <- c()
TimePoint <- c()
hpd_lower <- c()
hpd_upper <- c()

for (i in 1:length(log_files)) {
  medians <- medians_by_population[[i]]
  hpd_intervals <- hpd_by_population[[i]]
  
  for (j in 1:length(medians)) {
    Re_values <- c(Re_values, medians[j])
    hpd_lower <- c(hpd_lower, hpd_intervals[1, j])
    hpd_upper <- c(hpd_upper, hpd_intervals[2, j])
    Population <- c(Population, pop_names[i])
    TimePoint <- c(TimePoint, F_times[j])
    if (j %in% summer_indices) {
      Season <- c(Season, "Summer")
    } else {
      Season <- c(Season, "Winter")
    }
  }
}

# Create the data frame
data <- data.frame(
  Re = Re_values,
  Population = factor(Population),
  Season = factor(Season),
  TimePoint = TimePoint,
  hpd_lower = hpd_lower,
  hpd_upper = hpd_upper
)

# Calculate and print mean median Re values and standard error for each population and season
stats_by_pop_season <- data %>%
  group_by(Population, Season) %>%
  summarise(
    mean_median_Re = mean(Re, na.rm = TRUE),
    sd_median_Re = sd(Re, na.rm = TRUE),
    n = sum(!is.na(Re)),  # Count the number of non-NA values
    sd_median_Re = sd_median_Re / sqrt(n)  # Calculate standard error
  )


cat("Mean Median Re Values and Standard Deviations Across Time Points:\n")
print(stats_by_pop_season)

# Identify the specific winter periods where Re > 1
# For each population, find winter periods where median Re > 1
for (i in 1:length(log_files)) {
  pop_data <- data %>%
    filter(Population == pop_names[i], Season == "Winter")
  
  # Identify time points where median Re > 1
  over1_indices <- which(pop_data$Re > 1)
  
  if (length(over1_indices) > 0) {
    cat("\nPopulation:", pop_names[i], "has median Re > 1 during the following winter periods:\n")
    for (idx in over1_indices) {
      TimePoint_over1 <- pop_data$TimePoint[idx]
      Re_over1 <- pop_data$Re[idx]
      hpd_lower_over1 <- pop_data$hpd_lower[idx]
      hpd_upper_over1 <- pop_data$hpd_upper[idx]
      cat("Winter of", round(TimePoint_over1, 2), ":\n")
      cat("Median Re:", Re_over1, "\n")
      cat("95% HPD Interval: [", hpd_lower_over1, ",", hpd_upper_over1, "]\n")
    }
  } else {
    cat("\nPopulation:", pop_names[i], "did not have median Re > 1 during any winter periods.\n")
  }
}

