import argparse
import subprocess
import json
from tabulate import tabulate

def json_to_table(json_data):
    vResults = json_data['validationResults']
    drugResistance = json_data['drugResistance']
    dbResults = ""

    # Printing the validation results
    # print("Validation Results")
    for v in vResults:
        dbResults += "\n" + v['level'] + ": " + v['message'] + "\n"

    # Printing drug resistance results
    for d in drugResistance:
        dbResults += "\n\nDrug resistance interpretation: " + d['gene']['name']
        dbResults += "\nVersion: " + d['version']['text'] + "("  + d['version']['publishDate'] + ")"

        drugClass, drug, score, partialScores, text = [], [], [], [], []
        ctr = 0

        # Store the relevant information into a new dictionary
        for d2 in d['drugScores']:
            drugClass.append(d2['drugClass']['name'])
            drug.append(d2['drug']['displayAbbr'])
            score.append(d2['score'])
            partialScores.append(d2['partialScores'])
            text.append(d2['text'])
            ctr+=1

        results = {
            'drugClass' : drugClass,
            'drug' : drug,
            'score' : score,
            'text' : text
        }

        # Display Drug Resistance Results Table
        #display(pd.DataFrame(results))
        dbResults += "\n" + tabulate(results, headers='keys', tablefmt='psql') + "\n"

        # Printing partial scores, mutations, and comments
        for p in partialScores:
            if(p):
                mutations = p[0]['mutations']

                for m in mutations:
                    dbResults += f"\nMutation: {m['text']}\nType: {m['primaryType']}\nComments:\n"

                    for c in m['comments']:
                        dbResults += c['text'] + "\n"

    return dbResults

def create_sierra_cmd(aavf_file):
    cmd = "sierrapy mutations "

    for line in aavf_file:
        if line.startswith("#"):
            continue
        chrom, gene, pos, ref, alt, _, alt_freq, coverage, info = line.strip().split("\t")
        cmd = cmd + " " + gene + ":" + ref + pos + alt

    return cmd

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--aavf", required=True, help="AAVF input to query")
    parser.add_argument("--out", required=True, help="TXT output for rug resistance report")
    args = parser.parse_args()

    with open(args.aavf) as aavf_file, open(args.out, 'w') as txt_file:
        cmd = create_sierra_cmd(aavf_file)

        # Run the sierra command
        output = subprocess.check_output(cmd, shell=True)
        output_json_string = output.decode("utf-8")

        # Convert json output to more readable format
        output_json = json.loads(output_json_string)
        output_table = json_to_table(output_json)
        # print(output_table)
        txt_file.write(output_table)


if __name__ == '__main__':
    main()
