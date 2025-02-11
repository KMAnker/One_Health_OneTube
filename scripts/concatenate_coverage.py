import os
import pandas as pd
import re
import sys

# Global variable to store the user's choice for overwriting files
overwrite_existing_files = None

# Function to extract sample, method, OT, primers and segment from the Reference_Name column
def extract_sample_method_segment(reference_name):
    type_segment = re.search(r"(?<=_).*", reference_name).group(0)
    parts_segment = type_segment.split("_")
    sample_method = re.search(r"^[^_]+", reference_name).group(0)
    parts_sample_method = sample_method.split("-")
    
    # Assign values or NA for missing parts
    host = parts_sample_method[0]
    sample = parts_sample_method[0] + "-" + parts_sample_method[1]
    extraction = parts_sample_method[2] if len(parts_sample_method) >= 3 else None
    OT = parts_sample_method[3] if len(parts_sample_method) >= 4 else None
    primers = parts_sample_method[4] if len(parts_sample_method) >= 5 else None
    segment = parts_segment[1] if len(parts_segment) >= 2 else None
    segment2 = parts_segment[2] if len(parts_segment) >= 3 else parts_segment[1]
    
    return {
        "host": host,
        "sample": sample,
        "extraction": extraction,
        "OT": OT,
        "primers": primers,
        "segment": segment,
        "segment2": segment2
    }

# Function to modify the "Reference_Name" column and add new columns
def modify_and_rename_coverage_file(file_path):
    global overwrite_existing_files

    # Get the new file name (appending "_2" before the extension)
    new_file_name = file_path.replace("-coverage.txt", "-coverage_2.txt")

    # Check if the new file already exists
    if os.path.exists(new_file_name) and not overwrite_existing_files:
        print(f"Using existing file: {new_file_name}")
        return
    
    # Read the file
    data = pd.read_csv(file_path, delimiter='\t')

    # Extract the parent folder name from the file path
    folder_name = os.path.basename(os.path.dirname(os.path.dirname(file_path)))

    # Add the folder name to the "Reference_Name" column
    data['Reference_Name'] = folder_name + "_" + data['Reference_Name']

    # Extract and add new columns
    extracted_columns = data['Reference_Name'].apply(extract_sample_method_segment).apply(pd.Series)
    data = pd.concat([extracted_columns, data], axis=1)
    data = data.drop(columns=['Reference_Name'])
    
    data.to_csv(new_file_name, sep='\t', index=False)

# Function to recursively find all -coverage.txt files
def find_coverage_files(base_path):
    coverage_files = []
    for root, dirs, files in os.walk(base_path):
        # Extract the components of the current path
        path_parts = root.split(os.sep)
        # Check if the path has the desired structure and folder names
        if 'run' in path_parts:
            for file in files:
                if file.endswith("-coverage.txt") and not file.startswith("._"):
                    coverage_files.append(os.path.join(root, file))
    return coverage_files

# Function to check if any -coverage_2.txt files exist
def check_existing_files(coverage_files):
    global overwrite_existing_files

    # Check if any -coverage_2.txt files exist
    existing_files = [file.replace("-coverage.txt", "-coverage_2.txt") for file in coverage_files if os.path.exists(file.replace("-coverage.txt", "-coverage_2.txt"))]

    if existing_files:
        # Ask the user whether to overwrite existing files or use them
        response = input("Some -coverage_2.txt files already exist. Do you want to overwrite? (y/n): ")
        overwrite_existing_files = response.lower() == 'y'

# Function to find all -coverage.txt files, modify the "Reference_Name" column, and rename the files
def main(base_dir):
    coverage_files = find_coverage_files(base_dir)
    check_existing_files(coverage_files)
    if len(coverage_files) == 0:
        print("No -coverage.txt files found.")
    else:
        print("Found", len(coverage_files), "coverage.txt files.")
        for file_path in coverage_files:
            modify_and_rename_coverage_file(file_path)

    # Get the list of -coverage_2.txt files
    coverage_2_files = [file.replace("-coverage.txt", "-coverage_2.txt") for file in coverage_files]

    # Concatenate the modified files into coverage_all.txt
    with open(coverage_2_files[0], 'r') as first_file:
        header = first_file.readline()

    with open('coverage_all.txt', 'w') as concatenated_file:
        concatenated_file.write(header)
        for file_path in coverage_2_files:
            with open(file_path, 'r') as file:
                next(file)  # Skip the header
                for line in file:
                    concatenated_file.write(line)

    print("Concatenation complete. The result is saved in 'coverage_all.txt'")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python concatenate_coverage.py <base_dir>")
        sys.exit(1)

    main(sys.argv[1])