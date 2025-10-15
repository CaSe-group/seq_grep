#!/usr/bin/env python3
import os
import sys
import gzip
import re
import csv
import argparse
from typing import Dict, Set, List, Tuple
from collections import Counter

def search_barcodes(ids_file_path: str, fastq_dir_path: str, output_csv_path: str):
    """
    Searches FASTQ.gz files in a directory for specific read IDs and records the filename
    where each ID was found.

    Args:
        ids_file_path: Path to the .txt file containing read IDs (one per line).
        fastq_dir_path: Path to the directory containing .fastq.gz or .fq.gz files.
        output_csv_path: Path for the output CSV file (ID, Filename).
    """

    # --- REGEX DEFINITIONS (Header Parsing) ---
    # 1. Regex to extract the read ID from the FASTQ header (up to the first space/tab), without the leading '@'
    id_regex = re.compile(r'^@([^ \t]+)')
    # 2. Regex to handle potential read number suffixes (/1, /2) for ID normalization
    read_suffix_regex = re.compile(r'/[12]$')

    # 1. Load target IDs
    print(f"Reading target IDs from: {ids_file_path}")
    target_ids: Set[str] = set()
    try:
        with open(ids_file_path, 'r') as f:
            for line in f:
                read_id = line.strip()
                if read_id:
                    # Normalize the input ID (remove @ and /1 or /2 suffixes)
                    if read_id.startswith('@'):
                        read_id = read_id[1:]
                    read_id = read_suffix_regex.sub('', read_id)
                    target_ids.add(read_id)
    except FileNotFoundError:
        print(f"Error: IDs file not found at {ids_file_path}")
        sys.exit(1)

    if not target_ids:
        print("Warning: No IDs found in the input file. Exiting.")
        return

    # Initialize results dictionary. Keys are the normalized, non-prefixed IDs from the input.
    # The value will now store the filename.
    results: Dict[str, str] = {id_: 'NA' for id_ in target_ids}
    total_target_ids = len(target_ids)

    # Use a mutable set of IDs remaining to be found for efficiency.
    ids_to_find = target_ids.copy()

    # 2. Iterate through FASTQ files
    print(f"Searching FASTQ files in directory: {fastq_dir_path}")

    try:
        # Filter for compressed FASTQ files
        fastq_files = [f for f in os.listdir(fastq_dir_path) if f.endswith('.fastq.gz') or f.endswith('.fq.gz')]
    except FileNotFoundError:
        print(f"Error: FASTQ directory not found: {fastq_dir_path}")
        sys.exit(1)

    if not fastq_files:
        print(f"Warning: No .fastq.gz or .fq.gz files found in {fastq_dir_path}. Outputting 'NA' for all IDs.")

    for filename in fastq_files:
        # Stop searching if all IDs have been found
        if not ids_to_find:
            print("All target IDs found. Stopping search early.")
            break

        print(f"Processing: {filename}...")

        # The filename itself is the result to be recorded for any IDs found in this file.

        full_path = os.path.join(fastq_dir_path, filename)

        try:
            # Use 'rt' (read text) mode for gzip.open
            with gzip.open(full_path, 'rt', encoding='utf-8', errors='ignore') as fastq_file:

                # Process every 4th line (the header line)
                for line_count, line in enumerate(fastq_file, 1):
                    if line_count % 4 == 1 and line.startswith('@'):

                        # Extract the full ID from the FASTQ header
                        id_match = id_regex.search(line)
                        if not id_match:
                            continue

                        extracted_id = id_match.group(1).strip()

                        # Normalize the extracted ID (removes /1 or /2 suffix)
                        read_id = read_suffix_regex.sub('', extracted_id)

                        # Check if this base ID is one of our target IDs
                        if read_id in ids_to_find:
                            # Record the result using the filename
                            results[read_id] = filename
                            ids_to_find.remove(read_id)

        except Exception as e:
            print(f"An error occurred while reading {filename}: {e}")
            continue


    # 3. Write results to CSV
    print(f"\nWriting results to: {output_csv_path}")

    # Prepare the data list for CSV writer. Prepend '@' to the ID for the output.
    output_data: List[Tuple[str, str]] = [(f"@{id_}", results[id_]) for id_ in target_ids]

    try:
        with open(output_csv_path, 'w', newline='') as csvfile:
            csv_writer = csv.writer(csvfile)
            csv_writer.writerow(['ID', 'Filename']) # Updated header
            csv_writer.writerows(output_data)

        print(f"Successfully saved {len(output_data)} results.")

    except Exception as e:
        print(f"An error occurred while writing the CSV file: {e}")
        sys.exit(1)

    # 4. Summarize results and print to terminal
    print("\n--- Summary ---")

    # The summary now counts how many IDs were found in each filename.
    filename_counts = Counter(results.values())
    found_count = total_target_ids - filename_counts.get('NA', 0)
    percentage_found = (found_count / total_target_ids) * 100 if total_target_ids > 0 else 0
    
    # Sort the summary output by filename (alphabetically), handling 'NA' explicitly
    for filename, count in sorted(filename_counts.items()):
        if total_target_ids > 0:
            percentage_of_total = (count / total_target_ids) * 100
        else:
            percentage_of_total = 0.0

        if filename == 'NA':
            # Percentage for 'NA' is calculated relative to total IDs processed
            print(f"  {filename} (IDs Not Found): {count} ({percentage_of_total:.2f}%)")
        else:
            # Percentage for filename is calculated relative to total IDs processed
            print(f"  {filename}: {count} ({percentage_of_total:.2f}%)")
    # New section end

    print(f"-----------------------")
    print(f"  Total IDs Processed: {total_target_ids}")
    print(f"  IDs Found: {found_count}")
    print(f"  Percentage Found: {percentage_found:.2f}%")
    print(f"-----------------------")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Search for specific read IDs in FASTQ.gz files and record the filename where each ID was found.",
        formatter_class=argparse.RawTextHelpFormatter
    )

    # Positional Arguments
    parser.add_argument(
        'ids_file',
        type=str,
        help="Path to the .txt file containing target read IDs (one per line)."
    )

    parser.add_argument(
        'fastq_dir',
        type=str,
        help="Path to the directory containing .fastq.gz or .fq.gz files."
    )

    # Removed 'barcode_regex' argument

    # Optional Argument
    parser.add_argument(
        '-o', '--output',
        type=str,
        default='ID_to_filename_output.csv', # Changed default output name
        help="Path for the output CSV file (default: ID_to_filename_output.csv)."
    )

    # Parse arguments
    args = parser.parse_args()

    # Call the main function with parsed arguments
    search_barcodes(args.ids_file, args.fastq_dir, args.output)
