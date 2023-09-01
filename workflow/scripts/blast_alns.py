import argparse
import pandas as pd
import pysam

def get_overlap(qlist, slist):
    start1 = min(qlist)
    end1 = max(qlist)
    start2 = min(slist)
    end2 = max(slist)

    overlap = max(0, min(end1, end2) - max(start1, start2))
    #length = qend-qstart + send-sstart
    length = max(end1, end2) - min(start1, start2)
    # print("Overlap ", overlap)
    # print("Length ", length)
    # print("Return", overlap/length)

    return overlap*1.0/length

def is_same_pair(blast, index, read):
    return (blast[3].iloc[index]['qseqid'][-1] == '1' and read.flag < 128) or (blast[3].iloc[index]['qseqid'][-1] == '2' and read.flag >= 128)

def compare_alns(blast, aligner):
    print("in compare " + aligner[2])    
    sortedsam = aligner[3].replace(".sam", ".sorted.sam")
    pysam.sort("-n", "-o", sortedsam, aligner[3]) # sort samfile
    samfile = pysam.AlignmentFile(sortedsam, "r")

    correct = 0
    incorrect = 0
    total = 0
    mapped = 0
    unmapped = 0

    for read in samfile.fetch():
        if read.is_mapped:
            row = []
            key = ""
            if read.flag < 128: # is first pair
                key = read.query_name+"/1"
            else: # is second pair
                key = read.query_name+"/2"

            if key in blast[3]:
                mapped += 1
                row = blast[3][key]
                # print("found")
                # print(key, "\t", read.query_name, aligner[0], aligner[2])
            
                qstart1 = row['qstart']
                qend1 = row['qend']
                sstart1 = row['sstart']
                send1 = row['send']

                qstart2 = read.query_alignment_start
                qend2 = qstart2 + read.infer_query_length(always=False)
                sstart2 = read.reference_start
                send2 = read.reference_end

                # start, end
                qlist = [min(qstart1, qend1), max(qstart1, qend1)]
                slist = [min(qstart2, qend2), max(qstart2, qend2)]
                # print("Reference", qlist, slist)
                result_q = get_overlap(qlist, slist)

                qlist = [min(sstart1, send1), max(sstart1, send1)]
                slist = [min(sstart2, send2), max(sstart2, send2)]
                # print("Query", qlist, slist)
                result_s = get_overlap(qlist, slist)

                # print(result_q, result_s, "\n")
                if result_q >= 0.95 and result_s >= 0.95:
                    correct = correct + 1
                else:
                    incorrect = incorrect + 1

        else:
            unmapped += 1
            
        total+=1
        if total >= len(blast[3]):
            break
    print("Unmapped: ", unmapped)
    print("Mapped: ", mapped)
    print("Total: ", total)
    x = 1.0 * incorrect/mapped
    y = 1.0 * correct/mapped
    
    samfile.close()
    # print("Count", yes, correct, incorrect, total)
    return x, y

def read_blast(blast, file, sample, param, aligner):
    print("in read blast")
    data = pd.read_csv(file, sep='\t', names=["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"])
    data.sort_values('qseqid', inplace=True)
    data.set_index('qseqid', inplace=True)
    data = data[~data.index.duplicated(keep='first')]
    result_df = data.to_dict('index')
    blast.append([sample, param, aligner, result_df])

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--inFile', nargs="*")
    parser.add_argument('--outFile', required=True)
    args = parser.parse_args()

    blast = []
    last = []
    last_trained = []
    bowtie2 = []

    # separate input files according to the aligner used
    for file in args.inFile:
        my_str=file.split("/")
        country = my_str[2]
        sample = my_str[3]
        param = my_str[4]

        aligner = my_str[1].split("_")

        if aligner[1] == "trained":
            aligner = aligner[0] + " trained"
        else:
            aligner = aligner[0]

        if(aligner == "blast"):
            read_blast(blast, file, sample, param, aligner)
        elif(aligner == "last"):
            last.append([sample, param, aligner, file, country])
        elif(aligner == "last trained"):
            last_trained.append([sample, param, aligner, file, country])
        elif(aligner == "bowtie2"):
            bowtie2.append([sample, param, aligner, file, country])

    # x is incorrect, y is correct
    blast_bowtie_x = []
    blast_bowtie_y = []
    blast_last_x = []
    blast_last_y = []
    blast_lasttrained_x = []
    blast_lasttrained_y = []

    out = open(args.outFile, 'a')

    for i in range(len(blast)):
        # blast vs last
        x, y = compare_alns(blast[i], last[i])
        out.write(str(x) + '\t' + str(y) + '\t' + last[i][0] + '\t' + "LAST" + '\t' + last[i][4])
        out.write('\n')
        # blast_last_x.append(x)
        # blast_last_y.append(y)

        # blast vs last trained
        x, y = compare_alns(blast[i], last_trained[i])
        out.write(str(x) + '\t' + str(y) + '\t' + last_trained[i][0] + '\t' + "LAST (trained)" + '\t' + last_trained[i][4])
        out.write('\n')
        # blast_lasttrained_x.append(x)
        # blast_lasttrained_y.append(y)

        # blast vs bowtie2
        x, y = compare_alns(blast[i], bowtie2[i])
        out.write(str(x) + '\t' + str(y) + '\t' + bowtie2[i][0] + '\t' + "Bowtie2" + '\t' + bowtie2[i][4])
        out.write('\n')
        # blast_bowtie_x.append(x)
        # blast_bowtie_y.append(y)

    data = [blast_bowtie_x, blast_bowtie_y, blast_last_x, blast_last_y, blast_lasttrained_x, blast_lasttrained_y]

    # for d in data:
    #     for val in d:
    #         out.write(str(val) + '\t')
    #     out.write('\n')

    out.close()