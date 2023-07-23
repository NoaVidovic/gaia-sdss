import csv

# Input file path
input_file = "../data/stripe82calibStars_v4.2.dat"

# Output CSV file path
output_file = "../data/nonvariables.csv"

# Open the input file in read mode
with open(input_file, "r") as infile:
    # Read the lines of the file
    lines = infile.readlines()

# Extract the data and column names
data = []
for line in lines:
    if line.startswith("CALIBSTARS"):
        # Split the line by spaces
        columns = line.strip().split()
        # Extract the relevant data for each star (skip the first element, which is "CALIBSTARS_ID")
        data.append(columns[1:])

# Define the column names for the CSV file
column_names = [
    "ID", "RA", "Dec", "RArms", "Decrms", "Ntot", "Ar",
    "Nobs_u", "mmed_u", "mmu_u", "msig_u", "mrms_u", "mchi2_u",
    "Nobs_g", "mmed_g", "mmu_g", "msig_g", "mrms_g", "mchi2_g",
    "Nobs_r", "mmed_r", "mmu_r", "msig_r", "mrms_r", "mchi2_r",
    "Nobs_i", "mmed_i", "mmu_i", "msig_i", "mrms_i", "mchi2_i",
    "Nobs_z", "mmed_z", "mmu_z", "msig_z", "mrms_z", "mchi2_z"
]

# Write the data to the output CSV file with comma delimiters
with open(output_file, "w", newline="") as outfile:
    writer = csv.writer(outfile)
    # Write the header row with column names
    writer.writerow(column_names)
    # Write the data rows
    for idx, row in enumerate(data, start=1):
        writer.writerow([f"{idx}"] + row)

print(f"Data has been successfully saved to '{output_file}'.")
