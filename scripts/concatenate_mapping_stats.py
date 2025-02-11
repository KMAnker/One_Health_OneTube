import os
import pandas as pd
import sys

# Function to find mapping_stats.xls files with specific subdirectory structure
def find_mapping_stats_files(base_path):
    mapping_stats_files = []
    for root, dirs, files in os.walk(base_path):
        path_parts = root.split(os.sep)
        if 'run' in path_parts:
            for file in files:
                if file == "mapping_stats.xls":
                    mapping_stats_files.append(os.path.join(root, file))
    return mapping_stats_files

def main(base_dir):
    # Find all mapping_stats.xls files
    mapping_stats_files = find_mapping_stats_files(base_dir)

    # Create a new file with the header from the first file
    with open(mapping_stats_files[0], 'r') as first_file:
        header = first_file.readline()

    with open('mapping_stats_all.txt', 'w') as concatenated_file:
        concatenated_file.write(header)
        for file_path in mapping_stats_files:
            with open(file_path, 'r') as file:
                next(file)  # Skip the header
                for line in file:
                    concatenated_file.write(line)

    print("Concatenation complete. The result is saved in 'mapping_stats_all.txt'")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python concatenate_mapping_stats.py <base_dir>")
        sys.exit(1)

    main(sys.argv[1])