import os
import re
import pandas as pd
from collections import defaultdict
import sys

# Function to find all segments_validation.txt files
def find_segments_validation_files(base_path):
    segments_files = []
    for root, dirs, files in os.walk(base_path):
        path_parts = root.split(os.sep)
        if 'run' in path_parts:
            for file in files:
                if file == "segments_validation.txt":
                    segments_files.append(os.path.join(root, file))
    return segments_files


# Function to process a segments_validation.txt file and update the counts
def process_segments_validation_file(file_path, results):
    print(f"Processing file: {file_path}")
    with open(file_path, 'r') as file:
        for line in file:
            # Use regex to match the sample and method parts
            match = re.match(r'^([^-]+-[^-]+)-(.+?)_', line.strip())
            if match:
                sample = match.group(1)
                method = match.group(2)
                key = (sample, method)
                # Ensure the key is in the dictionary even if it has no 'Pass'
                if key not in results:
                    results[key] = 0
                if 'Pass' in line:
                    results[key] += 1


# Function to find all segments_validation.txt files
def main(base_dir):
    segments_files = find_segments_validation_files(base_dir)
    # Dictionary to hold the counts of passing segments
    results = defaultdict(int)
    # Process each segments_validation.txt file
    for file_path in segments_files:
        process_segments_validation_file(file_path, results)
    # Create a DataFrame from the results dictionary
    data = [{'Sample': sample, 'Method': method, 'passing_segments': count} for (sample, method), count in results.items()]
    df = pd.DataFrame(data)
    # Ensure the DataFrame has the expected columns before writing to file
    if not df.empty:
        df.to_csv('passing_segments_all.txt', sep='\t', index=False)
        print("Processing complete. The result is saved in 'passing_segments_all.txt'")
    else:
        print("No data to write. The DataFrame is empty.")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python get_passing_segments.py <base_dir>")
        sys.exit(1)

    main(sys.argv[1])
