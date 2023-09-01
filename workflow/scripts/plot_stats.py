import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse

parser = argparse.ArgumentParser(description='')
parser.add_argument('--tsvFile', nargs="*")
parser.add_argument('--outFile', required=True)
parser.add_argument('--countFile', required=True)
args = parser.parse_args()

last = pd.read_csv(args.tsvFile[0], sep='\t', header=0)
last_trained = pd.read_csv(args.tsvFile[1], sep='\t', header=0)
bowtie2 = pd.read_csv(args.tsvFile[2], sep='\t', header=0)

last["aligner"] = "LAST"
last_trained["aligner"] = "LAST (trained)"
bowtie2["aligner"] = "Bowtie2"

merged = pd.concat([last, last_trained, bowtie2])

merged.rename(columns = {'Sample ID':'Country'}, inplace = True)
merged['Parameters'] = merged['Parameters'].str.split("_", n=0, expand=True)[0]
merged.rename(columns = {'Parameters': 'Sample ID'}, inplace = True)
print(merged)

f = open(args.countFile, "w")
f.write(merged.to_csv(sep="\t", index=False))
f.close()

g = sns.barplot(data=merged, x='Sample ID', y='Alignment Percentage', hue='aligner')

# g.map(sns.barplot, 'sample_id', 'Alignment Percentage', 'aligner')
plt.xlabel('Sample ID')
plt.ylabel('Alignment Percentage')
plt.tight_layout()

# plt.show()

# plt.plot(last['Sample ID'], last['Alignment Percentage'], 'o',color='red', label='LAST')
# plt.plot(last_trained['Sample ID'], last_trained['Alignment Percentage'], 'o', color='orange', label='LAST (trained)')
# plt.plot(bowtie2['Sample ID'], bowtie2['Alignment Percentage'], 'o', color='blue', label='Bowtie2')

# plt.xlabel('Sample ID')
# plt.ylabel('Alignment Percentage')
# plt.legend(bbox_to_anchor=(1.05,1), loc='upper left')
# plt.tight_layout()

plt.savefig(args.outFile)

