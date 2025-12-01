#!/usr/bin/env python3
import gzip
import argparse
import os

def open_maybe_gzip(filename, mode="rt"):
    """
    Open a file that may be gzipped or plain text.
    """
    if filename.endswith(".gz"):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)

def make_fastq_unique(input_fastq, output_fastq):
    read_count = 1
    with open_maybe_gzip(input_fastq, "rt") as infile, open_maybe_gzip(output_fastq, "wt") as outfile:
        while True:
            header = infile.readline()
            if not header:
                break
            seq = infile.readline()
            plus = infile.readline()
            qual = infile.readline()

            # Split at first space to isolate read name
            header = header.strip()
            if " " in header:
                name, rest = header.split(" ", 1)
                name += f"_{read_count}"
                header = f"{name} {rest}\n"
            else:
                header = f"{header}_{read_count}\n"

            outfile.write(header)
            outfile.write(seq)
            outfile.write(plus)
            outfile.write(qual)

            read_count += 1

def main():
    parser = argparse.ArgumentParser(description="Make FASTQ read names unique")
    parser.add_argument("--input", "-i", required=True, nargs='+', help="Input FASTQ(.gz) file(s)")
    parser.add_argument("--sample_id", "-s", required=False, help="Sample ID (optional)")
    parser.add_argument("--outdir", "-o", required=True, help="Output directory")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # for idx, input_file in enumerate(args.input, start=1):
    #     if len(args.input) == 1:
    #         basename = f"{args.sample_id}_unique.fastq"
    #     else:
    #         basename = f"{args.sample_id}_{idx}_unique.fastq"

    #     output_file = os.path.join(args.outdir, basename)
    #     make_fastq_unique(input_file, output_file)

    basename = f"{args.sample_id}_unique.fastq"
    output_file = os.path.join(args.outdir, basename)
    make_fastq_unique(args.input[0], output_file)


if __name__ == "__main__":
    main()
