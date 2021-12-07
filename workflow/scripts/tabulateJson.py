from tabulate import tabulate
import json

if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--jsonFile', required=True)
    parser.add_argument('--outFile', required=True)
    args = parser.parse_args()

    data = json.load(open(args.jsonFile))
    vResults = data['validationResults']
    drugResistance = data['drugResistance']

    dbResults = ""

    # Printing the validation results
    # print("Validation Results")
    for v in vResults:
        dbResults += "\n" + v['level'] + ": " + v['message'] + "\n"

    # Printing drug resistance results
    for d in drugResistance:
        dbResults += "\n\nDrug resistance interpretation: " + d['gene']['name']
        dbResults += "\nVersion: " + d['version']['text'] + "("  + d['version']['publishDate'] + ")"

        drugClass = []
        drug = []
        score = []
        partialScores = []
        text = []
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
        dbResults += "\n"
        dbResults += tabulate(results, headers = 'keys', tablefmt = 'psql')
        dbResults += "\n"

        # Printing partial scores, mutations, and comments
        for p in partialScores:
            if(p):
                mutations = p[0]['mutations']

                for m in mutations:
                    dbResults += "\n"
                    dbResults += "Mutation: " + m['text'] + "\n"
                    dbResults += "Type: " + m['primaryType'] + "\n"
                    dbResults += "Comments:\n"

                    comments = m['comments']
                    for c in comments:
                        dbResults += c['text'] + "\n"

    #print(dbResults)

    f = open(args.outFile, 'w')
    f.write(dbResults)
    f.close()
    #print (output)
