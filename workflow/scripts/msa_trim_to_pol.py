import argparse
from Bio import SeqIO

def map_hxb2_coords_to_columns(hxb2_seq, start_pos, end_pos):
    alignment_start, alignment_end = None, None
    nt_count = 0
    for i, nt in enumerate(hxb2_seq):
        if nt != '-':
            nt_count += 1
        if nt_count == start_pos and alignment_start is None:
            alignment_start = i
        if nt_count == end_pos:
            alignment_end = i
            break
    return alignment_start, alignment_end

def trim_msa(input_file, output_file, hxb2_id, start, end, length_offset):
    start_pos = start + length_offset
    end_pos = end + length_offset

    # Read MSA
    records = list(SeqIO.parse(input_file, "fasta"))

    # Find HXB2 sequence
    hxb2_seq = None
    for r in records:
        if hxb2_id in r.id:
            hxb2_seq = r.seq
            break

    if not hxb2_seq:
        raise Exception(f"HXB2 ID '{hxb2_id}' not found in MSA")

    # Map to alignment columns
    alignment_start, alignment_end = map_hxb2_coords_to_columns(hxb2_seq, start_pos, end_pos)

    if alignment_start is None or alignment_end is None:
        raise Exception("Could not find full coordinate range in HXB2 sequence")

    print(f"Trimming alignment columns: {alignment_start} to {alignment_end} (HXB2 positions {start_pos}-{end_pos})")

    # Trim all sequences
    trimmed_records = []
    for r in records:
        r.seq = r.seq[alignment_start:alignment_end + 1]
        trimmed_records.append(r)

    # Write output
    SeqIO.write(trimmed_records, output_file, "fasta")
    print(f"Trimmed MSA saved to: {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Trim MSA to HXB2 pol region based on coordinates.")
    parser.add_argument("--input", required=True, help="Input MSA FASTA file")
    parser.add_argument("--output", required=True, help="Output trimmed FASTA file")
    parser.add_argument("--hxb2_id", default="K03455", help="Sequence ID for HXB2 in MSA (default: 'K03455')")
    parser.add_argument("--start", type=int, default=2253, help="Start position of region (default: 2253)")
    parser.add_argument("--end", type=int, default=5096, help="End position of region (default: 5096)")
    parser.add_argument("--offset", type=int, default=0, help="Optional offset to extend trimming range (default: 0)")

    args = parser.parse_args()
    trim_msa(args.input, args.output, args.hxb2_id, args.start, args.end, args.offset)
