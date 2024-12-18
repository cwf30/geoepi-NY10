from Bio import SeqIO
import csv
# FILEPATH: /Users/cwcf/Documents/ORISE/data/PMID26210076/process_PMID26210076.py

input_file = "PMID26210076_sequence.gb"
output_file = "PMID26210076_summary.csv"

# Open the input file and parse the records
records = SeqIO.parse(input_file, "genbank")

# Create a list to store the summarized features
summary = []

# Define the column names
column_names = ["UID", "accession", "collection_date", "isolation_source", "isolation_taxa", "city", "county", "state", "latitude", "longitude", "pseudo_lat", "pseudo_lat", "data_source"]

# Iterate over each record
for record in records:
    accession = record.id  # Extract the accession number
    accession = accession.split(".")[0]  # Remove the version number
    # Iterate over each feature in the record
    for feature in record.features:
        # Check if the feature is a source feature
        if feature.type == "source":
            # Extract the relevant information from the feature
            isolation_taxa = feature.qualifiers.get("host", [""])[0]
            latlong = feature.qualifiers.get("lat_lon", [""])[0]
            latitude = latlong.split(" ")[0]
            longitude = '-' + latlong.split(" ")[2]
            collection_date = feature.qualifiers.get("collection_date", [""])[0]
            # Append the summarized information to the list
            summary.append(["",accession, collection_date,"", isolation_taxa,"", "","",latitude, longitude,"","", 'PMID26210076'])

# Write the summarized features to the output CSV file
with open(output_file, "w", newline="") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(column_names)  # Write the header
    writer.writerows(summary)  # Write the data
