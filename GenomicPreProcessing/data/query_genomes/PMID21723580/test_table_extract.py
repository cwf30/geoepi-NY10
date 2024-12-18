import tabula

# Path to your file
file = "/Users/cwcf/Documents/ORISE/data/PMID21723580/NIHMS306066-supplement-01.pdf"

#convert pdf into csv
tabula.convert_into(file, "NIHMS306066-supplement-01.csv", output_format="csv", pages='all')
