#!/usr/bin/env python
# coding: utf-8



from __future__ import print_function
import vcf 
import sierrapy
import json
import requests
import subprocess
import os



aminoAcids ={   
#Phenylalanine    
"TTT": "F",
"TTC": "F",
#Leucine    
"TTA": "L",     
"TTG": "L",
"CTT": "L",
"CTC": "L",
"CTA": "L",
"CTG": "L", 
#Isoleucine
"ATT": "I",
"ATC": "I", 
"ATA": "I", 
#Methioinine    
"ATG": "M", 
#Valine  
"GTT": "V", 
"GTC": "V",
"GTA": "V",
"GTG": "V",
#Serine
"TCT": "S",
"TCC": "S",
"TCA": "S",
"TCG": "S",
"AGT": "S",
"AGC": "S",    
#Proline  
"CCT": "P",
"CCC": "P",
"CCA": "P",
"CCG": "P",
#Threonine
"ACT": "T",
"ACC": "T",
"ACA": "T",
"ACG": "T",      
#Alanine
"GCT": "A",
"GCC": "A",
"GCA": "A",
"GCG": "A",
#Tyrosine
"TAT": "Y",
"TAC": "Y",
#Histidine
"CAT": "H",
"CAC": "H",
#Glutamine
"CAA": "Q",
"CAG": "Q",
#Asparagine
"AAT": "N",
"AAC": "N",
#Lysine
"AAA": "K",
"AAG": "K",
#Aspartic acid
"GAT": "D",
"GAC": "D",
#Glutamic acid
"GAA": "E",    
"GAG": "E",
#Cysteine
"TGT": "C",
"TGC": "C",
#Tryptophan
"TGG": "W",
#Arginine
"CGT": "R",
"CGC": "R",
"CGA": "R",
"CGG": "R", 
"AGA": "R",
"AGG": "R",     
#Glycine
"GGT": "G",
"GGC": "G",
"GGA": "G",
"GGG": "G",    
#stop
"TAA": "*",
"TAG": "*",  
"TGA": "*"         
}




def getCharacter(position):
    fp = open('sequence.fasta', 'r')
    line_num = int(position//70 + 1) # add one to skip the first line
    char_num = int(position%70)
    print("Line Number: ", line_num, "Character Number: ", char_num)
    for i, line in enumerate(fp):
        if i == line_num:
            return line[char_num]
        



def getCodon(start,position,alt):
    adjusted = position - start
    print("Adjusted",adjusted)
    #start
    if adjusted%3 == 0:
        print("start")
        codon.append(alt)#codon.append(getCharacter(position))
        codon.append(getCharacter(position+1))
        codon.append(getCharacter(position+2))
    #middle
    elif adjusted%3 == 1:
        print("middle")
        codon.append(getCharacter(position-1))
        codon.append(alt)#codon.append(getCharacter(position))
        codon.append(getCharacter(position+1))
    #end
    elif adjusted%3 ==2:
        print("end")
        codon.append(getCharacter(position-2))
        codon.append(getCharacter(position-1))
        codon.append(alt)#codon.append(getCharacter(position))
    codon.append(adjusted+1) #position number in protease region. starts at 1
    return codon   



def getSequence(start,end):
    fp = open('sequence.fasta', 'r')
    print("Get Sequence")
    start_line = start//70 + 1
    end_line = end//70 + 1
    char_start = start%70
    char_end = end%70
    sequence = ""
    #print(start_line, end_line, char_start, char_end)
    for i, line in enumerate(fp):
        if i >= start_line and i <= end_line:
            if i == start_line:
                sequence = line[char_start:70]
            elif i == end_line:
                sequence = sequence + line[0:char_end+1]
            else:
                sequence = sequence + line[0:70]
            #print(sequence,char_end)
    fp.close()
    return sequence




def getGeneRegions():  
    
    #protease
    pStart = 2252
    pEnd = 2548
    protease = getSequence(pStart, pEnd)
    pf = open("protease.txt", "w")
    pf.write(protease)
    pf.close()
    #reverse transcriptase
    rtStart = 2549
    rtEnd = 4228 #3868 including rnase
    rTranscriptase = getSequence(rtStart, rtEnd)
    rf = open("reverse_transcriptase.txt", "w")
    rf.write(rTranscriptase)
    rf.close()
    iStart = 4229
    iEnd = 5092 #5095 stop codon is not included
    integrase = getSequence(iStart, iEnd)
    intf = open("integrase.txt", "w")
    intf.write(integrase)
    intf.close()
  



def getAminoAcids(inputF,outputF):

    inputFile = open(inputF, 'r')
    translation = ""
    sequence = ""
    outputFile= open(outputF, "w")
    for i, line in enumerate(inputFile):
        sequence = sequence + line
    seq_len = len(sequence)
    print("length of sequence",seq_len)
    for i in range(0,seq_len,3):
         #print("seq",sequence[i:i+3])
        translation = translation + aminoAcids.get(sequence[i:i+3])
    print("translated",translation)
    outputFile.write(translation)
    inputFile.close()
    outputFile.close()
    return translation



def getCombinations(base1,  base2, base3,reference):
# returns all possible combinations 
 
    mutations = []
    #print("Bases", base1, base2, base3)
    for i in range (len(base1)):
        for j in range (len(base2)):
            for k in range (len(base3)):
                codonStr = str(base1[i])+str(base2[j])+str(base3[k])
                #print("CODON:",codonStr)
                trans = aminoAcids.get(codonStr,"!")
                #print("TRANS:",trans, reference)
                if trans!="!" and trans!=reference: 
                    #add as a mutation only if it differs from the reference
                    if trans not in mutations:
                        #print("ALT codon",codonStr)
                        mutations.append(trans)
    return mutations        




def getMutations_mult(mut,refG): 
    #returns all possible amino acid mutations per position
    #protease
    mutList = []
    geneF = open(refG, 'r')
    ref_gene = [ch for ch in geneF.read()]
    codonct = 1
    #print("LENPROT",len(prot))
    for i in range(0,len(prot),3):
        #print("I",i,prot[i])
        #print("CODONCT ",codonct-1,ref_protease[codonct-1])
        #print("Base:",prot[i],prot[i+1],prot[i+2])
        combinations = getCombinations(mut[i],mut[i+1], mut[i+2],ref_gene[codonct-1])
        if len(combinations)>0:
            print("MUTATION:", codonct, ref_gene[codonct-1], "comb", combinations)
            mutList.append([codonct,ref_gene[codonct-1],combinations])
           
        codonct = codonct + 1
        
    print("M",mutList) 
    return mutList




def insertVariants_mult(vcfFile,pFileName, rtFileName, iFileName): 
    #gets all possible alternative bases per position
    vcf_reader = vcf.Reader(open(vcfFile,'r'))  
    pFile = open(pFileName, 'r')
    protease_seq = [ch for ch in pFile.read()]
    rtFile = open(rtFileName, 'r')
    rt_seq = [ch for ch in rtFile.read()]   
    iFile = open(iFileName, 'r')
    i_seq = [ch for ch in iFile.read()]
    #protease
    for record in vcf_reader:
       
        if record.POS >=2253 and record.POS <2550: #protease
            #print("variant",record.POS,record.POS-2253,record.REF,record.ALT)
            #print(record)
            alt = record.ALT
            
            temp = []
            #temp.append(record.REF)
            #temp.append(protease_seq[record.POS-2253])
            for i in range(len(alt[0])):
                base = str(alt[0])[i]
                if base not in temp:
                    temp.append(base)
            protease_seq[record.POS-2253] = temp 
            #print( "prot",protease_seq[record.POS-2253])
        elif record.POS >=2550 and record.POS < 4229: #3869: #reverse transcriptase
            #print("variant",record.POS,record.POS-2550,record.REF,record.ALT)
            alt = record.ALT
            temp = []
           #temp.append(record.REF)
            #temp.append(rt_seq[record.POS-2550])
            for i in range(len(alt[0])):
                base = str(alt[0])[i]
                if base not in temp:
                    temp.append(base)
            rt_seq[record.POS-2550] = temp 
            #print(rt_seq[record.POS-2550])
        elif record.POS >=4230 and record.POS <=5093: #5096 minus stop codon #integrase
            #print("variant",record.POS,record.POS-4230,record.REF,record.ALT)
            alt = record.ALT
            temp = []
            #temp.append(record.REF)
            #temp.append(i_seq[record.POS-4230])
            for i in range(len(alt[0])):
                base = str(alt[0])[i]
                if base not in temp:
                    temp.append(base)
            i_seq[record.POS-4230] = temp 
            #print(i_seq[record.POS-4230])
            
    return protease_seq, rt_seq, i_seq



def createHIVDBRequest(prot,rt,integrase, out):
    cmd = 'sierrapy mutations '
    #protease
    print("PROT1 ",len(prot))
    for i in prot:
        print("PROT: ", i)
        for j in i[2]:
            #print("SIERRA: " , i[0],i[1], j)
            mut = "PR:"+i[1]+str(i[0])+j
            print("SIERRA", mut)
            cmd = cmd + " "+mut
    #RT
    for i in rt:
        #print("PROT: ", i)
        for j in i[2]:
        #print("SIERRA: " , i[0],i[1], j)
            mut = "RT:"+i[1]+str(i[0])+j
            print("SIERRA", mut)
            cmd = cmd + " "+mut
    #RT
    for i in integrase:
        #print("PROT: ", i)
        for j in i[2]:
            #print("SIERRA: " , i[0],i[1], j)
            mut = "RT:"+i[1]+str(i[0])+j
            print("SIERRA", mut)
            cmd = cmd + " "+mut      
    #protStr
    cmd = cmd + " -o "+out
    #cmd = 'sierrapy mutations PR:L10I -o output_171.json'
    print(cmd)
    print (subprocess.check_output(cmd, shell=True))


"""

cmd = "sierrapy mutations PR:K64V -o output.json"
print (subprocess.check_output(cmd, shell=True))





vcfFile = '/Users/jenniferureta/Desktop/bioinfo/HIV/sequence/HIV2-171/171_mult.vcf'
prot,rt,i = insertVariants_mult(vcfFile,"protease.txt","reverse_transcriptase.txt","integrase.txt")
refp = "protease_amino.txt"
refrt = "rt_amino.txt"
refi = "integrase_amino.txt"
protease_mut = getMutations_mult(prot,refp)
print("PROTEASE MUT",protease_mut)
rt_mut = getMutations_mult(rt,refrt)
print("RT MUT",rt_mut)
int_mut = getMutations_mult(i,refi)
print("INT MUT",int_mut)
createHIVDBRequest(protease_mut, rt_mut,int_mut,"output_171.json")

"""


def insertVariants(vcfFile,pFileName, rtFileName, iFileName): 
    #one mutation per position and returns amino acids
    vcf_reader = vcf.Reader(open(vcfFile,'r'))
    
    pFile = open(pFileName, 'r')
    protease_seq = [ch for ch in pFile.read()]
    #print("PROT", protease_seq)
    rtFile = open(rtFileName, 'r')
    rt_seq = [ch for ch in rtFile.read()]   
    #print("RT", rt_seq)
    iFile = open(iFileName, 'r')
    i_seq = [ch for ch in iFile.read()]
    #print("I", i_seq)
    #lofreq starts at 1 but our indexing starts at 0 so always subtract 1 from the position given
    
    for record in vcf_reader:
       
        if record.POS >=2253 and record.POS <2550: #protease
            print("variant",record.POS,record.POS-2253,record.REF,record.ALT)
            alt = record.ALT[0]
            protease_seq[record.POS-2253] = alt
        elif record.POS >=2550 and record.POS < 4229: #3869: #reverse transcriptase
            print("variant",record.POS,record.POS-2550,record.REF,record.ALT)
            alt = record.ALT[0]
            rt_seq[record.POS-2550] = alt
        elif record.POS >=4230 and record.POS <=5093: #5096 minus stop codon #integrase
            print("variant",record.POS,record.POS-4230,record.REF,record.ALT)
            alt = record.ALT[0]
            i_seq[record.POS-4230] = alt
            
   # print("PROT VAR\n",protease_seq) 
    #print("RT VAR\n",rt_seq) 
    #print("I VAR\n",i_seq) 
    return protease_seq, rt_seq, i_seq





def getMutations(pseq, rtseq, iseq, refP, refRT, refI):
    #protease
    
    pmutations=[]
    rtmutations=[]
    imutations=[]
    plen = len(pseq)
    ptranslation = ""
    for i in range(0,plen,3): #translation to amino acid
        codonStr = str(pseq[i])+str(pseq[i+1])+str(pseq[i+2])
        #print("CODON",pseq[i],pseq[i+1],pseq[i+2], codonStr)
        ptranslation = ptranslation + aminoAcids.get(codonStr,"!")
    #output file format POSITION REF ALT
    pf = open(refP, 'r')
    ref_protease = [ch for ch in pf.read()]
    plen = len(ptranslation)

    for i in range(0,plen,1):
        if ref_protease[i] != ptranslation[i]:
            print("M", i+1,ref_protease[i],  ptranslation[i])
            pmutations.append([i+1,ref_protease[i],ptranslation[i]])
    #reverse transcriptase
    rtlen = len(rtseq)
    rttranslation = ""
    for i in range(0,rtlen,3):
        codonStr = str(rtseq[i])+str(rtseq[i+1])+str(rtseq[i+2])
        #print("CODON",rtseq[i],rtseq[i+1],rtseq[i+2], codonStr)
        rttranslation = rttranslation + aminoAcids.get(codonStr,"!")
       
    #output file format POSITION REF ALT
    rf = open(refRT, 'r')
    ref_RT = [ch for ch in rf.read()]
    rtlen = len(rttranslation)

    for i in range(0,rtlen,1):
        if ref_RT[i] != rttranslation[i]:
            print("M", i+1,ref_RT[i],  rttranslation[i])
            rtmutations.append([i+1,ref_RT[i],rttranslation[i]])
    #integrase
    ilen = len(iseq)
    itranslation = ""
    for i in range(0,ilen,3):
        codonStr = str(iseq[i])+str(iseq[i+1])+str(iseq[i+2])
        #print("CODON",iseq[i],iseq[i+1],iseq[i+2], codonStr)
        itranslation = itranslation + aminoAcids.get(codonStr,"!")
    #output file format POSITION REF ALT
    intf = open(refI, 'r')
    ref_int = [ch for ch in intf.read()]
    ilen = len(itranslation)

    for i in range(0,ilen,1):
        if ref_int[i] != itranslation[i]:
            print("M", i+1,ref_int[i],  itranslation[i])
            imutations.append([i+1,ref_int[i],itranslation[i]])
    print(pmutations)
    print(rtmutations)
    print(imutations)
    return pmutations,rtmutations,imutations


"""
#code for getting mutations
getGeneRegions()
getAminoAcids("protease.txt", "protease_amino.txt")
getAminoAcids("reverse_transcriptase.txt", "rt_amino.txt")
getAminoAcids("integrase.txt", "integrase_amino.txt")
vcfFile = '/Users/jenniferureta/Desktop/bioinfo/HIV/sequence/HIV2-171/171_sub.vcf'
pseq,rtseq,iseq = insertVariants(vcfFile,"protease.txt","reverse_transcriptase.txt","integrase.txt")
pm,rm,im= getMutations(pseq,rtseq,iseq,"protease_amino.txt", "rt_amino.txt", "integrase_amino.txt")
"""

#connectToHIVDB()

"""
codon = []

proteaseStart = 2252 #original protease start is 2253 for index starting at 1
rtStart = 2549 #original 2550
integrasetStart = 4229 #original 4230
vcf_reader = vcf.Reader(open('/Users/jenniferureta/Desktop/bioinfo/HIV/sequence/HIV2-171/171_sub.vcf','r'))
#lofreq starts at 1 but our indexing starts at 0 so always subtract 1 from the position given
mutations=[]
for record in vcf_reader:
    if record.POS >=2253 or record.POS <2550: #for protease gene 2253-2550
        mutation_position = int(record.POS)
        pCodon = getCodon(proteaseStart,mutation_position-1,record.ALT) 
        codonStr = pCodon[0]+pCodon[1]+pCodon[2]
        print(codonStr)
        print(aminoAcids.get(codonStr,0)) #0 if not found in amino acid table
#for reverse transcriptase gene 2550-3869


#for integrase gene 4230-5096

    print(record)
    print (record.POS)
"""




def convertSAMtoBAM(input_file, output_file):
    command = "samtools view -S -b "+input_file+" > "+output_file
    print(subprocess.call(command, shell = True))





def sortBAM(input_file, output_file):
    command = "samtools sort -o "+output_file+" " +input_file
    print(subprocess.call(command, shell = True))





def indexBAM(input_file):
    command = "samtools index "+ input_file
    print(subprocess.call(command, shell = True))




def performAlignment(file1, file2, out, log_file, index):
   
    command = "(bowtie2 --local -q -x "+ index + " -1 "+file1+ " -2 "+ file2 +" -S "+ out+" )2>> "+log_file
    print(command)
    subprocess.call(command, shell= True)




def callVariant(ref, inputf, outputf, log, threads):
    count =  str(threads)
    command = "(lofreq call-parallel --pp-threads +"+count+" -f"
    command = command + " " +ref +" -o " + outputf + " "+inputf+" )2>> "+log 
    subprocess.call(command, shell = True)




def buildIndex(inputf):
    command = "bowtie2-build "+inputf+ " hxb2"
    print(command)
    subprocess.call(command, shell= True)
    command = "samtools faidx "+inputf
    print(command)
    subprocess.call(command, shell= True)



in_dir = "/Users/jenniferureta/Desktop/bioinfo/HIV/sequence"

out_dir = "/Users/jenniferureta/Desktop/bioinfo/HIV/sequence_out"
logf = "/Users/jenniferureta/Desktop/bioinfo/HIV/sequence_out/file_test.log"
ref = "/Users/jenniferureta/Desktop/bioinfo/HIV/sequence.fasta"
#buildIndex(ref)
index = "hxb2"
for entry in os.scandir(in_dir):
    if entry.is_dir() and entry.name.startswith('HIV2'):
        file1 = in_dir+"/"+entry.name+"/"+entry.name +"_1.fastq"
        file2 = in_dir+"/"+entry.name+"/"+entry.name +"_2.fastq"
        out = out_dir+"/"+entry.name+".sam"
        filename = out_dir+"/"+entry.name
        print("OUT",out)
        #print(log_file)
        print("Path",os.path.exists(out_dir))
        if not os.path.exists(out_dir):
         # print("here")
            os.mkdir(out_dir,0o777)  
        #performAlignment(file1, file2, out, logf, index)
        outputf = out_dir+"/"+entry.name+".bam"
        #convertSAMtoBAM(out, outputf)
        #sortBAM(outputf, outputf)
        #indexBAM(outputf)
        vcf = out_dir+"/"+entry.name+".vcf"
        callVariant(ref, outputf, vcf, logf, 8)
        getGeneRegions()
        getAminoAcids("protease.txt", "protease_amino.txt")
        getAminoAcids("reverse_transcriptase.txt", "rt_amino.txt")
        getAminoAcids("integrase.txt", "integrase_amino.txt")
        prot,rt,i = insertVariants_mult(vcf,"protease.txt","reverse_transcriptase.txt","integrase.txt")
        refp = "protease_amino.txt"
        refrt = "rt_amino.txt"
        refi = "integrase_amino.txt"
        protease_mut = getMutations_mult(prot,refp)
        #print("PROTEASE MUT",protease_mut)
        rt_mut = getMutations_mult(rt,refrt)
        #print("RT MUT",rt_mut)
        int_mut = getMutations_mult(i,refi)
        #print("INT MUT",int_mut)
        outputf = out_dir+"/"+entry.name+".json"
        createHIVDBRequest(protease_mut, rt_mut,int_mut,outputf)



