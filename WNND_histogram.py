import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np

# =========================================
# 1. Data Loading and Preparation
# =========================================

# Load the WNV case count data
wnv_df = pd.read_csv('/path/to/human_neuroinvasive_wnv_2000-2021.csv')

# Load the population data with a specified encoding
pop_df = pd.read_csv('/path/to/co-est2019-alldata.csv', encoding='ISO-8859-1')

# Create the FIPS code by combining STATE and COUNTY columns
pop_df['fips'] = pop_df['STATE'].astype(str).str.zfill(2) + pop_df['COUNTY'].astype(str).str.zfill(3)

# Convert fips in wnv_df to string and ensure it's zero-padded
wnv_df['fips'] = wnv_df['fips'].astype(str).str.zfill(5)

# Filter relevant population columns (using only 2012 population estimates)
pop_df = pop_df[['fips', 'POPESTIMATE2012']]

# Merge the population data with the WNV case count data
wnv_df = pd.merge(wnv_df, pop_df, on='fips', how='left')

# Define the years of interest
years = [2011, 2012]

# Filter the data for the relevant years
wnv_filtered = wnv_df[wnv_df['year'].isin(years)]

# =========================================
# 2. Geospatial Processing
# =========================================

# Load your shapefile of US counties
shapefile_path = '/path/to/cb_2018_us_county_5m.shp'
gdf = gpd.read_file(shapefile_path)

# Ensure the FIPS codes are properly formatted (as strings and zero-padded if necessary)
gdf['GEOID'] = gdf['GEOID'].astype(str).str.zfill(5)
wnv_filtered['fips'] = wnv_filtered['fips'].apply(lambda x: str(x).zfill(5))

# Filter to keep only continental US (excluding Alaska, Hawaii, and other territories)
excluded_states = ['02', '15', '60', '66', '69', '72', '78']
continental_us = gdf[~gdf['STATEFP'].isin(excluded_states)]

# Merge the filtered shapefile with the filtered WNV case count data
gdf_merged = continental_us.merge(wnv_filtered, how='left', left_on='GEOID', right_on='fips')

# Re-project to a suitable projected CRS for accurate centroid calculation
gdf_merged = gdf_merged.to_crs(epsg=2163)

# Calculate the centroid coordinates in the projected CRS
gdf_merged['centroid_x'] = gdf_merged.geometry.centroid.x
gdf_merged['centroid_y'] = gdf_merged.geometry.centroid.y

# Convert the centroids back to the original geographic CRS
gdf_merged = gdf_merged.set_geometry(gpd.points_from_xy(gdf_merged['centroid_x'], gdf_merged['centroid_y'])).to_crs(epsg=4326)

# Extract the accurate longitude in familiar CRS
gdf_merged['longitude'] = gdf_merged.geometry.x

# Remove rows with NaN values in the 'count' column to avoid plotting errors
gdf_merged = gdf_merged.dropna(subset=['count', 'longitude'])

# =========================================
# 3. Filtering Longitudes Between -125 and -72
# =========================================

# Define longitude range
min_longitude = -123
max_longitude = -72

# Filter the data to include only longitudes between -125 and -72
gdf_filtered = gdf_merged[(gdf_merged['longitude'] >= min_longitude) & (gdf_merged['longitude'] <= max_longitude)]

# =========================================
# 4. Binning and Aggregation
# =========================================

# Define longitude bins (1-degree buckets) between -125 and -72
bins = np.arange(min_longitude, max_longitude + 1, 1)  # From -125 to -72 inclusive, step of 1

# Bin the longitudes
gdf_filtered['longitude_bin'] = pd.cut(gdf_filtered['longitude'], bins, include_lowest=True)

# Group by longitude_bin and year, then sum the cases
cases_by_longitude = gdf_filtered.groupby(['longitude_bin', 'year'])['count'].sum().unstack(fill_value=0)

# Extract the left edge (smallest longitude) of each bin for x-axis labels with °W notation
bin_labels = [f"{int(abs(bin.left))}°W" for bin in cases_by_longitude.index.categories]

# Update index to use the new bin labels based on the smallest longitude in each bin
cases_by_longitude.index = bin_labels


# =========================================
# 5. Plotting Side-by-Side Bars for Multiple Years
# =========================================

# Number of bins
n_bins = len(cases_by_longitude)

# Number of years
n_years = len(years)

# Positions on the x-axis for each bin
x = np.arange(n_bins)

# Width of each bar
bar_width = 0.8 / n_years  # Total width is 0.8, divided by number of years

# Create the plot
plt.figure(figsize=(20, 10))

# Define colors for each year (customize as needed)
colors = ['#444444', '#E30000']

# Plot bars for each year
for i, year in enumerate(years):
    plt.bar(x + (i - n_years/2) * bar_width + bar_width/2,
            cases_by_longitude[year],
            width=bar_width,
            color=colors[i],
            label=f'{year} Cases')

# Set x-axis ticks to align with bin midpoints and use the "°W" notation for labels
plt.xticks(ticks=x, labels=bin_labels, rotation=90, ha='center', fontsize=16)
plt.yticks(fontsize=16)
# Remove bounding box
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.gca().spines['left'].set_visible(False)
plt.gca().spines['bottom'].set_visible(False)

# Set labels and title with increased font size
plt.xlabel('Longitude (1° Bins, Larget Longitude Shown)', fontsize=24)
plt.ylabel('WNND Cases', fontsize=24)

# Add legend with doubled font size
plt.legend(title='Year', fontsize=18, title_fontsize=24)

# Adjust layout to prevent clipping
plt.tight_layout()
plt.savefig('/WNNDplots/WND_2012_histogram.png', dpi=300, bbox_inches='tight')

# Display the plot
plt.show()
