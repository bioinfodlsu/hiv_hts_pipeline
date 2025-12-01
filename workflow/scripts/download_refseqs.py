import os
import sys
from Bio import Entrez
from Bio import SeqIO
from pathlib import Path

# ACCESSIONS_PATH = Path(__file__).parent.parent.parent / "test_data" / "refseq_hiv1_m_accessions.txt"
ACCESSIONS_PATH = Path(__file__).parent.parent.parent / "test_data" / "m_subtypes_one_per.txt"
# ACCESSIONS_PATH = Path(__file__).parent.parent.parent / "test_data" / "simulated_skewed.txt"
# ACCESSIONS_PATH = Path(__file__).parent.parent.parent / "test_data" / "m_subtypes_3.txt"
DOWNLOADS_PATH = Path(__file__).parent.parent.parent / "test_data" / "refseqs"

def download_sequences(email):
    counter = 0
    with open(ACCESSIONS_PATH) as file:
        Entrez.email = email.strip()
        print(f"Using email: {Entrez.email}")
        print(f'Reading accessions from {ACCESSIONS_PATH}...')

        if not os.path.exists(ACCESSIONS_PATH):
            print(f"Error: Accessions file {ACCESSIONS_PATH} does not exist.", file=sys.stderr)
            sys.exit(1)

        for line in file:
            acc = line.strip()
            counter += 1
            print(f'[{counter}] Downloading {acc}...')
            handle = Entrez.efetch(db='nucleotide', id=acc, rettype='fasta', retmode='text')
            
            try:
                os.makedirs(DOWNLOADS_PATH)
            except FileExistsError:
                pass
            
            with open(f'{DOWNLOADS_PATH}/{acc}.fasta', 'w') as out_handle:
                out_handle.write(handle.read())
            handle.close()
            print(f'Downloaded {acc} successfully.')


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <email>", file=sys.stderr)
        sys.exit(1)

    download_sequences(sys.argv[1])
    print("All sequences downloaded successfully.")
