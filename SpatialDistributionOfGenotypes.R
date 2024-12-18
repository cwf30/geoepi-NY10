# Load required libraries
library(tidyverse)
library(maps)
library(viridis)
library(lubridate)
library(scatterpie)
library(sf)

# Read the data
metadata <- read.csv("/path/to/metadata.csv")
polymorphisms <- read.csv("/path/to/polymorphisms_dense.csv")

# Merge the data
data <- merge(metadata, polymorphisms, by="accession")

# Function to parse dates with multiple formats
parse_date <- function(x) {
  parsed <- dmy(x, quiet = TRUE)
  if (is.na(parsed) && nchar(x) == 4 && grepl("^\\d{4}$", x)) {
    parsed <- ymd(paste0(x, "-01-01"))  # Assume January 1st for year-only dates
  }
  return(parsed)
}

# Convert collection_date to proper Date format
data$collection_date_parsed <- as.Date(sapply(data$collection_date, parse_date), origin = "1970-01-01")

# Filter for years 2012-2021
data_filtered <- data %>%
  filter(year(collection_date_parsed) >= 2012 & year(collection_date_parsed) <= 2021)

# Create genotype variable
data_filtered <- data_filtered %>%
  mutate(
    genotype = case_when(
      `NS2A.R188` == 'K' & `NS4B.I240` == 'M' ~ "Both mutations",
      `NS2A.R188` == 'K' & `NS4B.I240` != 'M' ~ "NS2A.R188 K only",
      `NS2A.R188` != 'K' & `NS4B.I240` == 'M' ~ "NS4B.I240 M only",
      `NS2A.R188` != 'K' & `NS4B.I240` != 'M' ~ "Neither mutation",
      TRUE ~ NA_character_
    )
  )

# Calculate counts of each genotype per location
location_summary <- data_filtered %>%
  group_by(latitude, longitude) %>%
  summarize(
    samples = n(),
    both_mutations = sum(genotype == "Both mutations", na.rm = TRUE),
    ns2a_only = sum(genotype == "NS2A.R188 K only", na.rm = TRUE),
    ns4b_only = sum(genotype == "NS4B.I240 M only", na.rm = TRUE),
    neither_mutation = sum(genotype == "Neither mutation", na.rm = TRUE),
    .groups = "drop"
  )

# Convert map data to sf object
usa_sf <- st_as_sf(maps::map("state", plot = FALSE, fill = TRUE))

# Define the CRS (Albers Equal-Area)
crs_proj <- "+proj=aea +lat_1=39 +lat_2=45 +lon_0=-96 +datum=WGS84 +units=km +no_defs"

# Transform the USA map to the desired CRS
states_sf <- st_transform(usa_sf, crs = crs_proj)

# Convert location data to sf object and transform
locations_sf <- st_as_sf(location_summary, coords = c("longitude", "latitude"), crs = 4326)
locations_sf <- st_transform(locations_sf, crs = crs_proj)

# Extract coordinates for plotting
location_summary$x <- st_coordinates(locations_sf)[, 1]
location_summary$y <- st_coordinates(locations_sf)[, 2]


# Create color palette with corrected names
palette <- c(
  "both_mutations" = "#7c1158",
  "ns2a_only" = "#00B7C7",
  "ns4b_only" = "#FFCF25",
  "neither_mutation" = "#73af48"
)

# Calculate radius for pie charts based on samples and create a separate size variable for legend
max_radius <- 200  # Adjust for visible differences in radius
location_summary <- location_summary %>%
  mutate(
    radius = sqrt(samples) / sqrt(max(samples)) * max_radius,
    size_legend = samples  # Separate variable to map to size for legend
  )
# Extract coordinates for plotting after transforming
location_summary$x <- st_coordinates(locations_sf)[, 1]
location_summary$y <- st_coordinates(locations_sf)[, 2]

# Adjust legend_data to better align with scatterpie sizes by trial and error
legend_data <- data.frame(
  x = 0,  # Dummy x value for positioning in legend
  y = 2500,  # Dummy y value for positioning in legend
  radius = c(5,20,40,500),  # Trial values that visually match the pie chart sizes
  size_label = c("Small", "Medium", "Large", "Very Large")  # Descriptive labels
)

# Create the map
ggplot() +
  # Add the USA map with light grey fill and dark grey borders
  geom_sf(data = usa_sf, fill = "#fcf9f7", color = "darkgrey") +
  # Add state boundaries in a slightly darker shade
  geom_sf(data = states_sf, fill = NA, color = "darkgrey", size = 0.3) +
  # Scatter pie charts for mutation data with explicit radius mapping
  geom_scatterpie(
    aes(x = x, y = y, r = radius),
    data = location_summary,
    cols = c("both_mutations", "ns2a_only", "ns4b_only", "neither_mutation"),
    color = NA, alpha = 0.8
  ) +
  # Invisible points for size legend only
  geom_point(
    data = legend_data,
    aes(x = x, y = y, size = radius),
    color = NA, fill = NA, shape = 21, alpha = 0.5
  ) +
  # Size scale adjusted to better match scatterpie radii
  scale_size_continuous(
    name = "Sample Size",
    limits = c(0.5,500),  # Manually adjust to match the scatterpie scale
    breaks = c(10,20,300,500),
    labels = c("Small", "Medium", "Large", "Very Large"),
    guide = guide_legend(
      override.aes = list(shape = 21, fill = "grey"),
      title = "Sample Size",
      label.position = "right"
    )
  ) +
  # Fill scale for the pie chart colors
  scale_fill_manual(
    values = palette,
    name = "Genotype",
    labels = c(
      "both_mutations" = "NY10 (Both mutations)",
      "ns2a_only" = "NS2A-R188K only",
      "ns4b_only" = "NS4B-I240M only",
      "neither_mutation" = "Neither mutation"
    )
  ) +
  # Coordinate system and theme
  coord_sf(crs = crs_proj) +
  theme_void() +
  labs(
    title = ""
  ) +
  theme(
    text = element_text(size = 16),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 22, face = "bold", margin = margin(b = 20)),
    plot.margin = margin(t = 0, r = 20, b = 0, l = 0, unit = "pt")
  )

# Save the plot
ggsave("/Users/cwcf/Documents/ORISE/data/viz/figures/wnv_pie_map.png", width = 16, height = 12, dpi = 300)
