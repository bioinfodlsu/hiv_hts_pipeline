import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--tsvFile', required=True)
    parser.add_argument('--outFile', required=True)
    args = parser.parse_args()

    sns.set_theme()

    hiv = pd.read_csv(args.tsvFile, sep='\t', header=None, names=['Incorrect', 'Correct', 'Sample ID', 'Aligner', 'Country'])

    sns.relplot(
        data = hiv,
        x = 'Incorrect',
        y = 'Correct',
        hue = 'Aligner',
        style = 'Country'
    )

    plt.savefig(args.outFile)
    
    # data = []
    # file = open(args.tsvFile, "r")
    
    # for line in file:
    #     data.append([float(x) for x in line.strip().split('\t')])
    # file.close()

    # print(data)

    # plt.xlim(0)
    # plt.ylim(0)

    # plt.plot(data[0], data[1], 'o', color='blue', label='Bowtie2')
    # plt.plot(data[2], data[3], 'o', color='red', label='LAST')
    # plt.plot(data[4], data[5], 'o', color='orange', label='LAST (trained)')

    # plt.xlabel('Fraction Incorrectly Aligned')
    # plt.ylabel('Fraction Correctly Aligned')
    # plt.legend(bbox_to_anchor=(1.05,1), loc='upper left')
    # plt.tight_layout()

    # plt.savefig(args.outFile)