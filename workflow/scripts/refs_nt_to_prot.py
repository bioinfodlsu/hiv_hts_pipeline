#!/usr/bin/env python3
import os
import subprocess
from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(description="Extract gag and pol sequences from references using HXB2 coordinates.")
parser.add_argument("--hxb2", required=True, help="Path to HXB2 reference FASTA")
parser.add_argument("--refs", required=True, help="Path to multi-FASTA of your references")
parser.add_argument("--outdir", default="extracted_proteins", help="Output directory for translated protein sequences")
args = parser.parse_args()

# ---- USER INPUT ----
hxb2_fasta = args.hxb2   
ref_fasta = args.refs    
tempdir = os.path.join(os.path.dirname(args.outdir), "temp")
outdir = os.path.join(os.path.dirname(args.outdir), "extracted_cds")
protdir = args.outdir     # Output folder
# ---------------------

# hxb2_fasta = "/Users/gcoe/Documents/GitHub/pangenomics-hiv/test_data/references/hxb2.fasta"        # HXB2 GenBank file (download with efetch or from LANL)
# ref_fasta = "/Users/gcoe/Documents/GitHub/pangenomics-hiv/test_out/SouthAfrica/mafft/concat_refs.fasta"   # Your nucleotide reference FASTA (multi-fasta supported)
# tempdir = "temp"          # Temp folder for intermediate files
# outdir = "extracted_cds"   # Output folder
# protdir = "extracted_proteins"

os.makedirs(tempdir, exist_ok=True)
os.makedirs(outdir, exist_ok=True)
os.makedirs(protdir, exist_ok=True)

# 1) Load HXB2 reference
record = SeqIO.read(hxb2_fasta, "fasta")

# Known coordinates (HXB2 K03455.1, 1-based inclusive)
coords = {
    "gag": (790, 2292),
    "pol": (2253, 5096)
}

# 2) Extract gag/pol subsequences from HXB2
hxb2_queries = {}
for gene, (start, end) in coords.items():
    subseq = record.seq[start-1:end]  # convert to 0-based python indexing
    query_file = f"{tempdir}/HXB2_{gene}.fasta"
    with open(query_file, "w") as f:
        f.write(f">HXB2_{gene}\n{subseq}\n")
    hxb2_queries[gene] = query_file

# 3) Make BLAST database of your refs
subprocess.run(["makeblastdb", "-in", ref_fasta, "-dbtype", "nucl", "-out", f"{tempdir}/refsdb"], check=True)

# 4) BLAST gag/pol against refs and extract best hits
for gene, query in hxb2_queries.items():
    out_tsv = os.path.join(outdir, f"{gene}_vs_refs.tsv")
    cmd = [
        "blastn", "-query", query, "-db", f"{tempdir}/refsdb",
        "-outfmt", "6 qseqid sseqid qstart qend sstart send sstrand pident length evalue bitscore",
        "-max_target_seqs", "5", "-perc_identity", "70"
    ]
    with open(out_tsv, "w") as out:
        subprocess.run(cmd, stdout=out, check=True)

    # Parse BLAST hits: keep best (highest bitscore) per subject
    hits = {}
    with open(out_tsv) as f:
        for line in f:
            parts = line.strip().split("\t")
            sid, sstart, send, bitscore = parts[1], int(parts[4]), int(parts[5]), float(parts[10])
            if sid not in hits or bitscore > hits[sid][0]:
                hits[sid] = (bitscore, sstart, send)

    # Extract subsequences
    for rec in SeqIO.parse(ref_fasta, "fasta"):
        if rec.id in hits:
            _, start, end = hits[rec.id]
            if start > end:
                start, end = end, start
            subseq = rec.seq[start-1:end]  # 1-based → 0-based
            outfile = os.path.join(outdir, f"{rec.id}_{gene}.fasta")
            with open(outfile, "w") as f:
                f.write(f">{rec.id}_{gene}\n{subseq}\n")

print("Finished extracting gag and pol sequences.")

# Translate extracted sequences to protein
for gene in coords.keys():
    for fasta in os.listdir(outdir):
        if fasta.endswith(f"_{gene}.fasta") and not fasta.endswith("_prot.faa"):
            inpath = os.path.join(outdir, fasta)
            outpath = os.path.join(protdir, fasta.replace(".fasta", "_prot.faa"))
            with open(outpath, "w") as out_f:
                for rec in SeqIO.parse(inpath, "fasta"):
                    # Translate full sequence (ignoring internal stops)
                    prot_seq = rec.seq.translate(to_stop=False)
                    # Optional: remove trailing stops
                    prot_seq = prot_seq.rstrip("*")
                    out_f.write(f">{rec.id}_prot\n{prot_seq}\n")

print("Finished translating protein sequences.")