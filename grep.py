#!/usr/bin/python3
import sys
import re
import gzip
import argparse
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from itertools import islice

# --- CONSTANTS ---
# Define ANSI color codes for highlighting
COLOR_RED = "\033[91m"
COLOR_BLUE = "\033[94m"
COLOR_GREEN = "\033[92m"
COLOR_YELLOW = "\033[93m"
COLOR_RESET = "\033[0m"

# Number of reads (4 lines per read) to send to a worker process at once
BATCH_SIZE = 50000 
# Update progress bar every N batches
PROGRESS_UPDATE_FREQUENCY = 10 

# IUPAC Degeneracy Map (for converting sequence to regex pattern)
IUPAC_MAP = {
    'R': '[A|G]', 'Y': '[C|T]', 'S': '[G|C]', 'W': '[A|T]',
    'K': '[G|T]', 'M': '[A|C]', 'B': '[C|G|T]', 'D': '[A|G|T]',
    'H': '[A|C|T]', 'V': '[A|C|G]', 'N': '[A|C|G|T]',
    'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T',
}

# IUPAC Complement Map (for generating reverse complement sequence)
IUPAC_COMPLEMENT = {
    'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
    'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W',
    'K': 'M', 'M': 'K', 'B': 'V', 'D': 'H',
    'H': 'D', 'V': 'B', 'N': 'N',
}

# --- Utility Functions ---

def get_complement_base(base: str) -> str:
    """Returns the complementary base for a given IUPAC code."""
    return IUPAC_COMPLEMENT.get(base.upper(), base)

def reverse_complement_sequence(sequence: str) -> str:
    """Calculates the reverse complement of a DNA sequence, handling IUPAC codes."""
    complement = (get_complement_base(base) for base in sequence)
    return "".join(reversed(list(complement)))


def convert_to_regex(sequence: str) -> str:
    """
    Converts a sequence string containing IUPAC ambiguity codes into a
    case-insensitive regular expression pattern, respecting degeneracy.
    """
    if not sequence:
        return ""
    pattern = [IUPAC_MAP.get(base.upper(), re.escape(base)) for base in sequence]
    return "".join(pattern)


def highlight_matches(line: str, required_matches: list[dict]) -> str:
    """
    Highlights all occurrences that match the required regex patterns,
    using the specific color assigned to each logical term.
    
    NOTE: In the worker process, 'required_matches' must be reconstructed to include 
    color information necessary for highlighting.
    """
    highlighted_line = line
    for match_data in required_matches:
        color = match_data['color']
        for compiled_pattern in match_data['compiled_patterns']:
            if not compiled_pattern: continue

            def color_replacer(m, color=color):
                return f"{color}{m.group(0)}{COLOR_RESET}"

            highlighted_line = compiled_pattern.sub(
                color_replacer,
                highlighted_line
            )
    return highlighted_line.strip()

# --- Parallel Processing Worker Function (Batch-based) ---

def _search_read_batch(read_batch: list[tuple], required_matches_lite: list[dict]):
    """
    Worker function executed by the ProcessPoolExecutor.
    Processes a batch of FASTQ reads for matches.

    Args:
        read_batch: A list of tuples [(header_id, sequence_line), ...].
        required_matches_lite: Simplified match data (compiled patterns, term, color).

    Returns:
        A tuple: (found_matches, individual_hits_counts)
        - found_matches: List of (header_id, highlighted_output) for full matches.
        - individual_hits_counts: Dict mapping {term: count} of individual hits in this batch.
    """
    
    found_matches = []
    individual_hits_counts = {m['term']: 0 for m in required_matches_lite}

    # Reconstruct the full match data needed for highlighting (includes color)
    required_matches_full = [
        {
            'term': m['term'],
            'color': m['color'],
            'compiled_patterns': m['compiled_patterns'],
        } for m in required_matches_lite
    ]

    for header_id, sequence_line in read_batch:
        line_stripped = sequence_line.strip()
        all_logical_terms_present = True
        
        for match_data in required_matches_lite:
            term = match_data['term']
            is_individual_hit = False

            for p_compiled in match_data['compiled_patterns']:
                if p_compiled.search(line_stripped):
                    is_individual_hit = True
                    break

            if is_individual_hit:
                individual_hits_counts[term] += 1
            else:
                all_logical_terms_present = False

        if all_logical_terms_present:
            highlighted_output = highlight_matches(line_stripped, required_matches_full)
            found_matches.append((header_id, highlighted_output))

    return (found_matches, individual_hits_counts)

# --- Parallel Handler ---

def _read_fastq_batches(f_handle, batch_size):
    """
    Generator to yield batches of (header_id, sequence_line) tuples.
    """
    while True:
        current_batch = []
        for _ in range(batch_size):
            # Read 4 lines at a time (FASTQ block)
            lines = list(islice(f_handle, 4))
            if not lines:
                if current_batch:
                    yield current_batch
                return
            
            # Line 1: Header
            header = lines[0].strip()
            header_id = header.split(None, 1)[0] if header.startswith('@') else "Unknown"

            # Line 2: Sequence
            sequence = lines[1]
            
            current_batch.append((header_id, sequence))

        if current_batch:
            yield current_batch


def _process_reads_in_parallel(filepath: str, required_matches: list[dict], max_print_lines: int):
    """
    Handles reading the file sequentially and distributing sequence lines
    to a worker pool for parallel processing. (For uncompressed files only).
    """
    MAX_PRINT_LINES = max_print_lines

    # Prepare a lightweight version of required_matches for pickling/passing to workers
    # We include 'color' here so the worker can perform highlighting if a full match occurs.
    required_matches_lite = [
        {
            'term': m['term'],
            'color': m['color'],
            'compiled_patterns': m['compiled_patterns'],
        } for m in required_matches
    ]

    opener = open
    mode = 'rt'

    print(f"--- Processing file: {filepath} (Parallel Search) ---")
    print(f"--- Note: Output of matching reads is limited to the first {MAX_PRINT_LINES} examples. ---")

    num_cores = max(2, os.cpu_count() or 1)
    print(f"Using {num_cores} worker processes with a batch size of {BATCH_SIZE} reads.")

    futures = []
    found_count = 0
    reads_processed = 0
    batch_counter = 0

    try:
        with opener(filepath, mode, encoding='utf-8') as f:
            # 1. Submit batches of read blocks to the process pool
            with ProcessPoolExecutor(max_workers=num_cores) as executor:
                for read_batch in _read_fastq_batches(f, BATCH_SIZE):
                    reads_processed += len(read_batch)
                    batch_counter += 1
                    futures.append(executor.submit(
                        _search_read_batch,
                        read_batch,
                        required_matches_lite
                    ))
                
                # 2. Collect results as they complete (as_completed)
                reads_completed_approx = 0
                for future in as_completed(futures):
                    try:
                        found_matches_in_batch, individual_hits_counts = future.result()

                        # Aggregate individual hits in the main process
                        for term, count in individual_hits_counts.items():
                            for match_data in required_matches:
                                if match_data['term'] == term:
                                    match_data['total_individual_count'] += count
                                    break
                        
                        # Handle full match printing (Printing is done sequentially in the main thread)
                        for header_id, highlighted_output in found_matches_in_batch:
                            found_count += 1
                            if found_count <= MAX_PRINT_LINES:
                                sys.stdout.write(f"\r[{header_id}]: {highlighted_output}\n")
                                sys.stdout.flush()
                            elif found_count == MAX_PRINT_LINES + 1:
                                sys.stdout.write(f"\r[... Output limited after {MAX_PRINT_LINES} reads to accelerate processing ...]\n")
                                sys.stdout.flush()
                        
                        # Update progress bar
                        reads_completed_approx += BATCH_SIZE # Approximate completion status
                        if batch_counter % PROGRESS_UPDATE_FREQUENCY == 0:
                            # Use total submitted reads for the final count
                            sys.stdout.write(f"\rReads processed: {reads_completed_approx:,} (Total submitted: {reads_processed:,})")
                            sys.stdout.flush()


                    except Exception as e:
                        print(f"\nError in worker process: {e}", file=sys.stderr)

        # Print final status
        sys.stdout.write(f"\rReads processed: {reads_processed:,} (100% complete)\n")
        return reads_processed, found_count

    except FileNotFoundError:
        print(f"\nError: The file '{filepath}' was not found.")
        return 0, 0
    except Exception as e:
        print(f"\nAn unexpected error occurred during parallel processing: {e}", file=sys.stderr)
        return 0, 0


def _grep_and_highlight_sequential(filepath: str, required_matches: list[dict], max_print_lines: int):
    """
    Searches a file sequentially (used for gzipped files).
    """
    MAX_PRINT_LINES = max_print_lines
    PROGRESS_UPDATE_FREQUENCY = 10000

    is_compressed = filepath.lower().endswith(('.gz', '.tgz'))
    opener = gzip.open if is_compressed else open
    mode = 'rt'

    try:
        with opener(filepath, mode, encoding='utf-8') as f:
            line_number = 0
            found_count = 0
            last_header_id = None
            
            print(f"--- Processing file: {filepath} (Sequential Search) ---")
            print(f"--- Note: Output of matching reads is limited to the first {MAX_PRINT_LINES} examples. ---")
            
            for line in f:
                line_number += 1
                line_stripped = line.strip()

                # --- FASTQ Header Logic: Capture the read ID (Line 1: @) ---
                if line_number % 4 == 1 and line.startswith('@'):
                    try:
                        last_header_id = line_stripped.split(None, 1)[0]
                    except IndexError:
                        last_header_id = "@Unknown"
                
                # --- Sequence Line Logic (Line 2: The actual read sequence) ---
                if line_number % 4 == 2:
                    reads_processed = line_number // 4
                    
                    # --- PROGRESS BAR UPDATE ---
                    if reads_processed % PROGRESS_UPDATE_FREQUENCY == 0:
                        sys.stdout.write(f"\rReads processed: {reads_processed:,}")
                        sys.stdout.flush()
                    # --- END PROGRESS BAR UPDATE ---

                    all_logical_terms_present = True
                    
                    for match_data in required_matches:
                        is_individual_hit = False
                        
                        for p_compiled in match_data['compiled_patterns']:
                            if p_compiled.search(line_stripped):
                                is_individual_hit = True
                                break
                        
                        if is_individual_hit:
                            match_data['total_individual_count'] += 1
                        else:
                            all_logical_terms_present = False
                    
                    if all_logical_terms_present:
                        found_count += 1
                        
                        if found_count <= MAX_PRINT_LINES:
                            identifier = last_header_id if last_header_id is not None else f"Line {line_number}"
                            highlighted_output = highlight_matches(line_stripped, required_matches)
                            
                            sys.stdout.write(f"\r[{identifier}]: {highlighted_output}\n")
                            sys.stdout.flush()
                        elif found_count == MAX_PRINT_LINES + 1:
                            sys.stdout.write(f"\r[... Output limited after {MAX_PRINT_LINES} reads to accelerate processing ...]\n")
                            sys.stdout.flush()

            # Final update
            final_reads_processed = line_number // 4
            sys.stdout.write(f"\rReads processed: {final_reads_processed:,} (100% complete)\n")
            return final_reads_processed, found_count

    except FileNotFoundError:
        print(f"\nError: The file '{filepath}' was not found.")
        return 0, 0
    except Exception as e:
        print(f"\nAn unexpected error occurred: {e}", file=sys.stderr)
        return 0, 0


def grep_and_highlight(filepath: str, required_matches: list[dict], max_print_lines: int):
    """
    Router function to select parallel or sequential processing.
    """
    if not required_matches:
        print("Error: Please provide at least one search term.")
        return

    # Determine if the file is compressed
    is_compressed = filepath.lower().endswith(('.gz', '.tgz'))

    if is_compressed:
        # GZIP files cannot be easily parallelized, use sequential streaming
        total_reads, found_count = _grep_and_highlight_sequential(filepath, required_matches, max_print_lines)
    else:
        # Use parallel processing for CPU-intensive regex matching
        total_reads, found_count = _process_reads_in_parallel(filepath, required_matches, max_print_lines)

    # --- FINAL SUMMARY ---
    print(f"\n--- Search Complete ---")
    
    # Calculate the percentage of reads that matched all required terms
    percentage_match = 0.0
    if total_reads > 0:
        percentage_match = (found_count / total_reads) * 100
        
    lines_printed = min(found_count, max_print_lines) 
    
    print(f"Total reads processed: {total_reads}")
    print(f"Total matching reads (ALL terms present, printed: {lines_printed}): {found_count} ({percentage_match:.2f}%)")

    print("--- Individual Term Hit Summary: ---")
    
    for match_data in required_matches:
        term = match_data['term']
        total_individual_count = match_data['total_individual_count']
        color = match_data['color']
        
        colored_total_count = f"{color}{total_individual_count}{COLOR_RESET}"
        colored_term = f"{color}{term}{COLOR_RESET}"
        
        print(f"  {colored_total_count} total hits for term: {colored_term}")
    print(f"-----------------------")
    # --- END SUMMARY ---


# --- Setup and Execution ---

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Grep-like tool to find lines containing ALL specified search terms in a file (including .gz). Supports IUPAC degeneracy and colors results. Optimized for speed using parallel processing for uncompressed files."
    )
    
    parser.add_argument(
        'filepath',
        type=str,
        help="Path to the input file (.txt, .fastq, or compressed .gz file)."
    )
    
    # Argument for controlling output limit
    parser.add_argument(
        '-l', '--limit',
        type=int,
        default=20, 
        help="Maximum number of matching reads to print to the terminal output (Default: 20)."
    )
    
    parser.add_argument(
        '-rc', '--reverse_complement',
        action='store_true',
        help="Also search for the reverse complement of every provided sequence."
    )
    
    color_groups = [
        ('red', COLOR_RED), 
        ('blue', COLOR_BLUE), 
        ('green', COLOR_GREEN), 
        ('yellow', COLOR_YELLOW)
    ]
    
    # Add 36 optional search terms (4 colors * 9 indices)
    args_to_collect = {}
    for color_name, _ in color_groups:
        for i in range(1, 10):
            flag_long = f'--{color_name}{i}'
            flag_short = f'-{color_name[0]}{i}'
            
            parser.add_argument(
                flag_short, flag_long, 
                type=str, 
                default='', 
                help=f"Required term {i} to be highlighted in {color_name.capitalize()}."
            )
            args_to_collect[flag_long.lstrip('--')] = color_name
    
    args = parser.parse_args()
    
    # 5. Collect, filter, convert terms, and handle Reverse Complement
    required_matches = []
    
    for arg_name, color_name in args_to_collect.items():
        term = getattr(args, arg_name)
        
        if term and term.strip():
            term_fwd = term.strip()
            
            color_code = next(c for n, c in color_groups if n == color_name)
            
            pattern_strings = [convert_to_regex(term_fwd)]
            
            if args.reverse_complement:
                term_rc = reverse_complement_sequence(term_fwd)
                pattern_strings.append(convert_to_regex(term_rc))
            
            # Compile all pattern strings once
            compiled_patterns = [
                re.compile(p, re.IGNORECASE) for p in pattern_strings
            ]

            required_matches.append({
                'term': term_fwd,
                'color': color_code,
                'compiled_patterns': compiled_patterns,
                'total_individual_count': 0
            })
    
    if not required_matches:
        parser.error("At least one search term (e.g., -r1, --blue2) must be provided.")
    
    # 6. Call the main function, passing the limit argument
    grep_and_highlight(args.filepath, required_matches, args.limit)
