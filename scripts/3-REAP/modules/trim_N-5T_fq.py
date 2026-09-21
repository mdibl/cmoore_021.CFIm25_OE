############ allow gz input, thanks for Dinghai Zheng and Luyang Wang. # 2025/3/25

import argparse
import os
import gzip
import re
from Bio.Seq import Seq

# Define the Fastq class
class Fastq:
    """ Represents a Fastq record """
    def __init__(self, name, seq, qual):
        self.name = name
        self.seq = seq
        self.qual = qual

    def get_name(self):
        return self.name

    def get_seq(self):
        return self.seq

    def get_qual(self):
        return self.qual

# Generator for reading Fastq files
def reader_fastq(infile):
    """Generates Fastq objects from a given file, supports .gz"""
    
    i = 0
    name, seq, qual = None, None, None
    open_func = gzip.open if infile.endswith('.gz') else open
    with open_func(infile, 'rt') as f:
        for line in f:
            i += 1
            curr_line = line.strip()
            if i % 4 == 1:
                name = curr_line[1:]
            elif i % 4 == 2:
                seq = curr_line
            elif i % 4 == 0:
                qual = curr_line
                yield Fastq(name, seq, qual)

# Function to trim random nucleotides and trailing Ts
def fastq_trim_Ts(infile, outfile, random_NT_len=6, min_t_length=8, max_mismatch=2):
    """
    Trim 5' T sequences with allowed mismatches.
    Supports gzipped input/output.
    """
    print(f"Processing {infile}...", flush=True)
    open_out = gzip.open if outfile.endswith('.gz') else open
    with open_out(outfile, 'wt') as outhandle:
        count = 0  # Count of reads written to output
        for fq in reader_fastq(infile):
            seq, qual = fq.get_seq(), fq.get_qual()

            # Trim random nucleotides at the start
            seq, qual = seq[random_NT_len:], qual[random_NT_len:]

            # Detect and trim T-rich regions with allowed mismatches
            t_region_end, mismatch_count = 0, 0
            for i, base in enumerate(seq):
                if base == 'T':
                    t_region_end = i + 1
                else:
                    mismatch_count += 1
                    if mismatch_count > max_mismatch:
                        break

            # Check if the T-rich region is long enough
            if t_region_end >= min_t_length:
                seq, qual = seq[t_region_end:], qual[t_region_end:]

            # Write trimmed sequence if it's long enough
            if len(seq) >= 18:
                count += 1
                outhandle.write(f"@{fq.get_name()}\n{seq}\n+\n{qual}\n")

    print(f"Processed {count} reads from {infile}. Output saved to {outfile}.", flush=True)

# Main function
def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Trim 5' random nucleotides and trailing Ts from a Fastq file.")
    parser.add_argument("--rawfastq", required=True, metavar='<rawfastq>',
                        help="Path to the input Fastq file (supports .fastq, .fq, .fastq.gz, .fq.gz).")
    parser.add_argument("--out_dir", required=True, metavar='<out_dir>',
                        help="Working out directory.")
    parser.add_argument("--random_NT_len", type=int, default=6, metavar='<random_NT_len>',
                        help="Length of random nucleotides to trim (default: 6).")

    # Parse arguments
    args = parser.parse_args()

    # Set parameters
    rawfastq, out_dir, random_NT_len = args.rawfastq, args.out_dir, args.random_NT_len

    # Create result directory if not exists
    os.makedirs(out_dir, exist_ok=True)

    # Change to the out directory
    os.chdir(out_dir)

    # Determine output file name
    output_file = os.path.join(out_dir, os.path.basename(rawfastq).replace('.fastq', '.5Ttrimmed.fastq').replace('.fq', '.5Ttrimmed.fq'))

    # Process Fastq file
    fastq_trim_Ts(infile=rawfastq, outfile=output_file, random_NT_len=random_NT_len)

    print(f"Processed file saved as: {output_file}", flush=True)

if __name__ == "__main__":
    main()