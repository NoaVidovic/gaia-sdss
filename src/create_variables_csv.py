import csv

input_file = "../data/stripe82candidateVar_v1.1.dat.gz"
output_file = "../data/variables.csv"

# Provided header from the data
header = [
    "ID", "ra", "dec", "P", "r", "ug", "gr", "ri", "iz",
    "gN", "gAmpl", "rN", "rAmpl", "iN", "iAmpl", "zQSO", "MiQSO"
]

# Function to parse the data and convert it to a list of lists
def parse_data(file_path):
    data = []
    with open(file_path, "r") as file:
        for line in file:
            # Skip lines starting with '#' and any leading/trailing whitespaces
            if not line.startswith("#"):
                line_data = line.strip().split()
                # Append non-empty lines to the data list
                if line_data:
                    data.append(line_data)
    return data


# Function to export the data to a CSV file
def export_to_csv(data, file_path):
    with open(file_path, "w", newline="") as file:
        csv_writer = csv.writer(file)
        # Write the header first
        csv_writer.writerow(header)
        # Write the data rows
        for row in data:
            csv_writer.writerow(row)


# Parse the data from the input file
parsed_data = parse_data(input_file)

# Export the parsed data to a CSV file with comma delimiter
export_to_csv(parsed_data, output_file)

print(f"The data has been exported to {output_file} successfully.")
