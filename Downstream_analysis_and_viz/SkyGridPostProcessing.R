library(ape)
library(phylodyn)
library(ggplot2)

# Define colors for the populations
colors <- c("#7C1158", "#11A579")  # Colors for the two lineages

# Load your first MCC tree
NY10_tree <- read.nexus("/path/to/NY10_SkyGrid.tree")

# Load your second MCC tree
WN02_tree <- read.nexus("/path/to/WN02_SkyGrid.tree")  

# Run BNPR_PS for each tree
BNPR_PS_results_NY10 <- BNPR_PS(data = NY10_tree, lengthout = 11)
BNPR_PS_resultsWN02 <- BNPR_PS(data = WN02_tree, lengthout = 16)

# Extract the relevant data for the first tree
plot_data_NY10 <- data.frame(
  time = BNPR_PS_results_NY10$x,
  effpopmean = BNPR_PS_results_NY10$effpopmean,
  effpop025 = BNPR_PS_results_NY10$effpop025,
  effpop975 = BNPR_PS_results_NY10$effpop975
)

# Extract the relevant data for the second tree
plot_data_WN02 <- data.frame(
  time = BNPR_PS_resultsWN02$x,
  effpopmean = BNPR_PS_resultsWN02$effpopmean,
  effpop025 = BNPR_PS_resultsWN02$effpop025,
  effpop975 = BNPR_PS_resultsWN02$effpop975
)

# Convert time to actual years for both datasets, based on approximate last sample dates (Oct 7 for NY10 and September 27 for WN02)
plot_data_NY10$year <- 2018.76 - plot_data_NY10$time
plot_data_WN02$year <- 2018.74 - plot_data_WN02$time

# Apply natural log transformation to effective population size values for both datasets
plot_data_NY10$log_effpopmean <- log(plot_data_NY10$effpopmean)
plot_data_NY10$log_effpop025 <- log(plot_data_NY10$effpop025)
plot_data_NY10$log_effpop975 <- log(plot_data_NY10$effpop975)

plot_data_WN02$log_effpopmean <- log(plot_data_WN02$effpopmean)
plot_data_WN02$log_effpop025 <- log(plot_data_WN02$effpop025)
plot_data_WN02$log_effpop975 <- log(plot_data_WN02$effpop975)

# Add a group identifier to each dataset
plot_data_NY10$group <- "Group II"
plot_data_WN02$group <- "WN02"

# Combine the two datasets
plot_data <- rbind(plot_data_NY10, plot_data_WN02)

# Plot with custom legend labels
ggplot() +
  # HPD interval as shaded area for each group
  geom_ribbon(data = plot_data, aes(x = year, ymin = log_effpop025, ymax = log_effpop975, fill = group), alpha = 0.2) +
  # Median line for each group
  geom_path(data = plot_data, aes(x = year, y = log_effpopmean, color = group), size = 1.2) +
  # Points for medians of each group
  geom_point(data = plot_data, aes(x = year, y = log_effpopmean, color = group), size = 2) +
  
  # Styling
  labs(x = "Year", y = expression(N[e] ~ "(Log Effective population size)"), color = "Group", fill = "Group") +
  theme_minimal(base_size = 15) +
  theme(
    axis.text.x = element_text(size = 12, color = "black"),  # Set x-axis text color to black
    axis.title.x = element_text(color = "black"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    axis.ticks = element_line(color = "black"),
    panel.grid = element_blank(),  # Remove background grid
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    legend.position = c(0.85, 0.30),  # Place legend inside plot area, top-right corner
    legend.background = element_rect(fill = "white", color="white")  # Add white background and border to legend
  ) +
  scale_x_continuous(breaks = c(2010, 2012, 2014, 2016, 2018), limits = c(2008, 2019)) +  # Restrict x-axis to 2009-2018
  geom_hline(yintercept = 0, color = "grey") +  # Baseline for natural log scale
  # Custom colors and legend labels
  scale_color_manual(values = colors, labels = c("Group II", "WN02")) +
  scale_fill_manual(values = colors, labels = c("Group II", "WN02"))