import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as mcolors
from matplotlib.colors import TwoSlopeNorm, LinearSegmentedColormap
import imageio
import os


# Set Helvetica font and larger default text size
plt.rcParams['font.family'] = 'Helvetica'
plt.rcParams['font.size'] = 24  # General text size
plt.rcParams['axes.titlesize'] = 36  # Title font size
plt.rcParams['axes.labelsize'] = 24  # Axis label font size
plt.rcParams['legend.fontsize'] = 18  # Legend font size

# Load the data
df = pd.read_csv('/path/to/human_neuroinvasive_wnv_2000-2021.csv')

# Load the population data
pop_df = pd.read_csv('/path/to/co-est2019-alldata.csv', encoding='ISO-8859-1')

# Create the FIPS code by combining STATE and COUNTY columns in the population data
pop_df['fips'] = pop_df['STATE'].astype(str).str.zfill(2) + pop_df['COUNTY'].astype(str).str.zfill(3)

# Convert fips in WNV data to string and ensure it's zero-padded
df['fips'] = df['fips'].astype(str).str.zfill(5)

# Merge the population data with the WNV case count data
df = pd.merge(df, pop_df, on='fips', how='left')

# Load your shapefile of US counties
gdf = gpd.read_file('/path/to/cb_2018_us_county_5m.shp')

# Ensure the FIPS codes are properly formatted (as strings and zero-padded if necessary)
gdf['GEOID'] = gdf['GEOID'].astype(str)

# Filter to keep only continental US (excluding Alaska, Hawaii, and other territories)
continental_us = gdf[~gdf['STATEFP'].isin(['02', '15', '60', '66', '69', '72', '78'])]

# Define the Albers Equal-Area projection
crs_proj = "+proj=aea +lat_1=39 +lat_2=45 +lon_0=-96 +datum=WGS84 +units=km +no_defs"

# Create a directory for saving individual year maps if it doesn't exist
output_dir = '/WNNDplots'
os.makedirs(output_dir, exist_ok=True)

# Initialize a list to store image file paths for GIF creation
image_files = []

# Initialize a previous year's incidence DataFrame
prev_year_incidence = None

# Loop through each year and calculate the change in incidence
for year in range(2001, 2020+1):
    # Filter the data for the specific year and create a copy
    df_year = df[df['year'] == year].copy()
    
    # Calculate incidence for the year using the 2012 population
    df_year['incidence'] = (df_year['count'] / df_year['POPESTIMATE2012']) * 100000
    
    if prev_year_incidence is not None:
        # Calculate the change in incidence from the previous year
        df_year = pd.merge(df_year, prev_year_incidence, on='fips', how='left', suffixes=('', '_prev'))
        df_year['incidence_change'] = df_year['incidence'] - df_year['incidence_prev']
    else:
        # No change for the first year, so just use incidence
        df_year['incidence_change'] = df_year['incidence']

    # Update the previous year's incidence
    prev_year_incidence = df_year[['fips', 'incidence']].copy()

    if year == 2001:
        continue  # Skip saving a map for 2001, since it's just used to calculate the change for 2002 onward

    # Merge the shapefile with the change in incidence data
    gdf_filtered = continental_us.merge(df_year, how='left', left_on='GEOID', right_on='fips')

    # Apply the Albers Equal-Area projection to the GeoDataFrame
    gdf_filtered = gdf_filtered.to_crs(crs_proj)

    # Create a custom colormap that transitions through white at 0
    colors = ["blue", "white", "red"]
    custom_cmap = LinearSegmentedColormap.from_list("custom_coolwarm", colors, N=256)

    # Continue with your existing code...

    # Plotting and color normalization
    fig, ax = plt.subplots(1, 1, figsize=(15, 10))

    # Calculate lower and upper caps, ensuring they are correctly ordered
    lower_cap = gdf_filtered['incidence_change'].quantile(0.01)  # Cap at the 1st percentile
    upper_cap = gdf_filtered['incidence_change'].quantile(0.99)  # Cap at the 99th percentile

    # Find the maximum absolute cap to make the color scale symmetrical
    max_abs = max(abs(lower_cap), abs(upper_cap))

    # Set symmetrical vmin and vmax
    vmin = -max_abs
    vmax = max_abs
    vcenter = 0

    norm = TwoSlopeNorm(vmin=vmin, vcenter=vcenter, vmax=vmax)

    # Plot with the custom colormap
    gdf_filtered.plot(column='incidence_change', ax=ax, cmap=custom_cmap, legend=False, norm=norm,
                    missing_kwds={"color": "white"})

    # Overlay the county boundaries
    gdf_filtered.boundary.plot(ax=ax, linewidth=0.5, edgecolor='0.75')

    # Add the color bar with the custom colormap
    sm = plt.cm.ScalarMappable(cmap=custom_cmap, norm=norm)
    sm._A = []
    cbar = fig.colorbar(sm, ax=ax, orientation='horizontal', aspect=20, pad=0, shrink=0.5)
    cbar.set_label(f"Change in WNND Incidence in {year} (from previous year)", fontsize=18)

    # Increase title font size
    plt.axis('off')
    plt.tight_layout()

# Save or show the figure as desired

    year_map_path = os.path.join(output_dir, f'wnv_incidence_change_{year}.png')
    plt.savefig(year_map_path)
    plt.close()

    # Append the file path to the list of images for GIF creation
    image_files.append(year_map_path)

# Create a multi-panel figure for all years (excluding 2001)
fig, axes = plt.subplots(nrows=4, ncols=5, figsize=(20, 15))

# Flatten the axes array for easy iteration
axes = axes.flatten()

for i, year in enumerate(range(2002, 2020+1)):
    # Load the saved image for each year
    img = plt.imread(image_files[i])

    # Plot the image in the corresponding subplot
    axes[i].imshow(img)
    axes[i].set_title(f'{year}', fontsize=18)
    axes[i].axis('off')

plt.tight_layout()
plt.show()


