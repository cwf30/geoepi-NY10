from Bio import SeqIO
import csv

input_file = "PMID30457531_sequence.gb"
output_file = "PMID30457531_summary.csv"

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
            country = feature.qualifiers.get("country", [""])[0]
            # Extract the state and county names from the country field
            
            if ":" in country:
                list = country.split(": ")[1].split(", ")
                if len(list) == 2:
                    if 'Collins' in list[1]:
                        state, city, county = list[0], "Fort Collins", "Larimer"
                    else:
                        state, city, county = list[0], "", list[1].split(" ")[0]
                elif len(list) == 3:
                    state, city, county = list[0], list[1], list[2].split(" ")[0]
                else:
                    state, city, county = list[0], "", ""
            print(f'state: {state}, city: {city}, county: {county}')
            collection_date = feature.qualifiers.get("collection_date", [""])[0]
            # Append the summarized information to the list
            if city != "" or county != "":
                summary.append(["", accession, collection_date, "", isolation_taxa, city, county, state, "", "", "", "", 'PMID30457531'])

# Write the summarized features to the output CSV file
with open(output_file, "w", newline="") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(column_names)  # Write the header
    writer.writerows(summary)  # Write the data