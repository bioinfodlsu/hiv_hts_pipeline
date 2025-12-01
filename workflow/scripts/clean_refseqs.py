from Bio import SeqIO
from pathlib import Path

# ACCESSIONS_PATH = Path(__file__).parent.parent.parent / "test_data" / "refseq_hiv1_m_accessions.txt"
ACCESSIONS_PATH = Path(__file__).parent.parent.parent / "test_data" / "m_subtypes_one_per.txt"
# ACCESSIONS_PATH = Path(__file__).parent.parent.parent / "test_data" / "m_subtypes_3.txt"
DOWNLOADS_PATH = Path(__file__).parent.parent.parent / "test_data" / "refseqs"
CLEANED_PATH = Path(__file__).parent.parent.parent / "test_data" / "refseqs_cleaned"

def clean_fasta(in_file, out_file):
    print(f"Cleaning FASTA from {in_file} to {out_file}...")
    allowed = set("ACGTNacgtn-")
    with open(out_file, "w") as out_handle:
        for record in SeqIO.parse(in_file, "fasta"):
            clean_seq = ''.join([b if b in allowed else 'N' for b in str(record.seq)])
            record.seq = clean_seq
            SeqIO.write(record, out_handle, "fasta")

if __name__ == "__main__":

    if not CLEANED_PATH.exists():
        CLEANED_PATH.mkdir(parents=True)

    for fasta_file in Path(DOWNLOADS_PATH).glob("*.fasta"):
        clean_fasta(fasta_file, Path(CLEANED_PATH) / fasta_file.name)

    print("Cleaning complete.")