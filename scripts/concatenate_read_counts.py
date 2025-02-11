import os
import pandas as pd
import sys

# Global variable to store the user's choice for overwriting files
overwrite_existing_files = None

# Function to modify the "Record" column and save as a new READ_COUNTS_2.txt file
def modify_and_rename_read_counts_file(file_path):
    global overwrite_existing_files

    # Get the new file name (appending "_2" before the extension)
    new_file_name = file_path.replace("READ_COUNTS.txt", "READ_COUNTS_2.txt")

     # Check if the new file already exists
    if os.path.exists(new_file_name) and not overwrite_existing_files:
        print(f"Using existing file: {new_file_name}")
        return

    # Read the file
    data = pd.read_csv(file_path, delimiter='\t')
    
    # Extract the parent folder name from the file path
    folder_name = os.path.basename(os.path.dirname(os.path.dirname(file_path)))
    
    # Add the folder name to the "Record" column
    data['Record'] = folder_name + "_" + data['Record']
    
    # Save the modified data to the new file
    data.to_csv(new_file_name, sep='\t', index=False)

# Function to find READ_COUNTS.txt files with specific subdirectory structure
def find_read_counts_files(base_path):
    read_counts_files = []
    for root, dirs, files in os.walk(base_path):
        path_parts = root.split(os.sep)
        if 'run' in path_parts:
            for file in files:
                if file == "READ_COUNTS.txt":
                    read_counts_files.append(os.path.join(root, file))
    return read_counts_files

# Function to check if any READ_COUNTS_2.txt files exist
def check_existing_files(read_counts_files):
    global overwrite_existing_files

    # Check if any READ_COUNTS_2.txt files exist
    existing_files = [file.replace("READ_COUNTS.txt", "READ_COUNTS_2.txt") for file in read_counts_files if os.path.exists(file.replace("READ_COUNTS.txt", "READ_COUNTS_2.txt"))]

    if existing_files:
        # Ask the user whether to overwrite existing files or use them
        response = input("Some READ_COUNTS_2.txt files already exist. Do you want to overwrite? (y/n): ")
        overwrite_existing_files = response.lower() == 'y'

# Function to find all READ_COUNTS.txt files, modify the "Record" column, and rename the files
def main(base_dir):
    read_counts_files = find_read_counts_files(base_dir)
    check_existing_files(read_counts_files)
    for file_path in read_counts_files:
        modify_and_rename_read_counts_file(file_path)

    # Get the list of READ_COUNTS_2.txt files
    read_counts_2_files = [file.replace("READ_COUNTS.txt", "READ_COUNTS_2.txt") for file in read_counts_files]

    # Create a new file with the header from the first file
    with open(read_counts_2_files[0], 'r') as first_file:
        header = first_file.readline()

    with open('READ_COUNTS_all.txt', 'w') as concatenated_file:
        concatenated_file.write(header)
        for file_path in read_counts_2_files:
            with open(file_path, 'r') as file:
                next(file)  # Skip the header
                for line in file:
                    concatenated_file.write(line)

    print("Concatenation complete. The result is saved in 'READ_COUNTS_all.txt'")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python concatenate_read_counts.py <base_dir>")
        sys.exit(1)

    main(sys.argv[1])