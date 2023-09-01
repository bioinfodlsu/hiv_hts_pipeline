import argparse
import subprocess

parser = argparse.ArgumentParser(description='')
parser.add_argument('--inFile', nargs="*")
parser.add_argument('--outFile', required=True)
args = parser.parse_args()

out = open(args.outFile, "w")

out.write("Sample ID\tParameters\tTotal Mapped\tTotal Unmapped\tTotal Reads\tAlignment Percentage\n")

results = []
reads = {}
for f in args.inFile:

    # cmd = "samtools flagstat " + f
    # cmd_output = subprocess.check_output(cmd, shell=True)
    # out.write(cmd_output + "\n\n")
    
    names = f.split('/')
    sample_id = names[2]
    params = names[3]

    cmd = "samtools view -c -F 4 " + f
    cmd_output = subprocess.check_output(cmd, shell=True)
    mapped_count = cmd_output.decode("utf-8").strip()

    if "bowtie2" not in f:
        file = f.replace("last_trained", "bowtie2") if "last_trained" in f else f.replace("last", "bowtie2")
        # print(f)
        # print(file)
        cmd = "samtools view -c " + file
    else:
        cmd = "samtools view -c " + f
    
    cmd_output = subprocess.check_output(cmd, shell=True)
    alns_count = cmd_output.decode("utf-8").strip()

    unmapped_count = str(int(alns_count) - int(mapped_count))

    aln_pct = '%.3f'%(float(mapped_count)/float(alns_count)*100)

    out.write(sample_id + "\t" + 
              params + "\t" + 
              mapped_count + "\t" + 
              unmapped_count + "\t" +
              alns_count + "\t" + 
              aln_pct + "\n")

out.close()
