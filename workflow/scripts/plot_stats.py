import pandas as pd
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser(description='')
parser.add_argument('--tsvFile', nargs="*")
parser.add_argument('--outFile', required=True)
args = parser.parse_args()

last = pd.read_csv(args.tsvFile[0], sep='\t', header=0)
last_trained = pd.read_csv(args.tsvFile[1], sep='\t', header=0)
bowtie2 = pd.read_csv(args.tsvFile[2], sep='\t', header=0)

plt.plot(last['Sample ID'], last['Alignment Percentage'], 'o',color='red', label='LAST')
plt.plot(last_trained['Sample ID'], last_trained['Alignment Percentage'], 'o', color='orange', label='LAST (trained)')
plt.plot(bowtie2['Sample ID'], bowtie2['Alignment Percentage'], 'o', color='blue', label='Bowtie2')

plt.xlabel('Sample ID')
plt.ylabel('Alignment Percentage')
plt.legend(bbox_to_anchor=(1.05,1), loc='upper left')
plt.tight_layout()

plt.savefig(args.outFile)

