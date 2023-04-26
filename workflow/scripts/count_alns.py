import argparse
import subprocess

parser = argparse.ArgumentParser(description='')
parser.add_argument('--inFile', nargs="*")
parser.add_argument('--outFile', required=True)
args = parser.parse_args()

out = open(args.outFile, "w")

out.write("Sample ID\tParameters\tTotal Mapped\tTotal Unmapped\tTotal Reads\tAlignment Percentage\n")

results = []
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

    cmd = "samtools view -c -f 4 " + f
    cmd_output = subprocess.check_output(cmd, shell=True)
    unmapped_count = cmd_output.decode("utf-8").strip() 

    cmd = "samtools view -c " + f
    cmd_output = subprocess.check_output(cmd, shell=True)
    alns_count = cmd_output.decode("utf-8").strip()

    aln_pct = '%.3f'%(float(mapped_count)/float(alns_count)*100)

    out.write(sample_id + "\t" + 
              params + "\t" + 
              mapped_count + "\t" + 
              unmapped_count + "\t" +
              alns_count + "\t" + 
              aln_pct + "\n")

out.close()
