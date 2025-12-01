import pysam
import argparse

HXB2_REF = "K03455.1"
HXB2_FASTA = "/Users/gcoe/Documents/GitHub/pangenomics-hiv/test_data/references/hxb2.fasta"

def read_fasta_to_dict(fasta_path, hxb2_path=HXB2_FASTA):
    """Read FASTA file(s) and return a dictionary of {header: sequence}."""
    haplotypes = {}

    def _parse_fasta(path):
        current, seq_chunks = None, []
        with open(path, 'r') as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                if line.startswith('>'):
                    if current is not None:
                        haplotypes[current] = ''.join(seq_chunks)
                    current = line[1:].split()[0]
                    seq_chunks = []
                else:
                    seq_chunks.append(line)
            if current is not None:
                haplotypes[current] = ''.join(seq_chunks)

    _parse_fasta(fasta_path)
    _parse_fasta(hxb2_path)

    return haplotypes


def build_match_mismatch_cigar(qseq, ref_seq):
    """Generate CIGAR tuples for matches (=), mismatches (X), insertions (I), deletions (D)."""
    min_len = min(len(qseq), len(ref_seq))
    if min_len == 0:
        return []

    cig_runs = []
    cur_op = 7 if qseq[0] == ref_seq[0] else 8
    cur_len = 1

    for i in range(1, min_len):
        op = 7 if qseq[i] == ref_seq[i] else 8
        if op == cur_op:
            cur_len += 1
        else:
            cig_runs.append((cur_op, cur_len))
            cur_op, cur_len = op, 1
    cig_runs.append((cur_op, cur_len))

    if len(qseq) > min_len:
        cig_runs.append((1, len(qseq) - min_len))  # insertion
    if len(ref_seq) > min_len:
        cig_runs.append((2, len(ref_seq) - min_len))  # deletion

    return cig_runs


def recalibrate_read(read, haplotypes, target_ref=HXB2_REF):
    """Reassign a read to the target reference and build a new CIGAR if sequences differ."""
    qseq = read.query_sequence or ''
    ref_name = read.reference_name

    if ref_name not in haplotypes:
        return None

    new_ref_seq = haplotypes[target_ref]
    ref_slice = new_ref_seq[read.reference_start:read.reference_start + len(qseq)]

    new_read = pysam.AlignedSegment()
    new_read.query_name = read.query_name
    new_read.query_sequence = qseq
    new_read.flag = read.flag & ~0x4  # ensure mapped
    new_read.reference_start = read.reference_start
    new_read.mapping_quality = read.mapping_quality
    new_read.query_qualities = read.query_qualities

    # Build new CIGAR if mismatch
    if ref_slice != qseq:
        new_read.cigartuples = build_match_mismatch_cigar(qseq, ref_slice)
    else:
        new_read.cigartuples = read.cigartuples

    return new_read


def recalibrate_sam(input_sam, output_sam, haplotypes, target_ref=HXB2_REF):
    """Process all reads from SAM and output recalibrated SAM."""
    with pysam.AlignmentFile(input_sam, "r") as sam_in:
        header = sam_in.header.to_dict()

        # Add HXB2 to header if missing
        if target_ref not in [sq['SN'] for sq in header.get('SQ', [])]:
            header['SQ'].append({'SN': target_ref, 'LN': len(haplotypes[target_ref])})

        with pysam.AlignmentFile(output_sam, "w", header=header) as sam_out:
            target_tid = sam_out.get_tid(target_ref)

            for read in sam_in.fetch():
                new_read = recalibrate_read(read, haplotypes, target_ref)
                if new_read:
                    # assign correct reference in the output header context
                    new_read.reference_id = target_tid
                    sam_out.write(new_read)

    print(f"✅ Recalibrated SAM written to: {output_sam}")


def main():
    parser = argparse.ArgumentParser(
        description="Reassign reads from a SAM file to a target reference (HXB2) using haplotype FASTA sequences."
    )
    parser.add_argument(
        "--haplotypes",
        nargs="+",
        required=True,
        help="Path(s) to FASTA file(s) containing reference haplotypes (can be multiple)."
    )
    parser.add_argument(
        "--input",
        required=True,
        help="Path to input SAM file."
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Path to output recalibrated SAM file."
    )

    args = parser.parse_args()

    # Load haplotypes from one or more FASTA files
    haplotypes = {}
    for fasta_path in args.haplotypes:
        haplotypes.update(read_fasta_to_dict(fasta_path))
    print(f"Loaded {len(haplotypes)} reference sequences from {len(args.haplotypes)} FASTA file(s)")

    recalibrate_sam(args.input, args.output, haplotypes)


if __name__ == "__main__":
    main()
