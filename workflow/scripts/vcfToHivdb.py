import os
import vcf
#from pysam import VariantFile
import subprocess
import datetime

#coordinates of the 3 genes in 1-based coordinates (because VCF is 1-based)
#protease
pStart = 2253
pEnd = 2549
#reverse-transcriptase
rtStart = 2550
rtEnd = 4229
#integrase
iStart = 4230
iEnd = 5093

#HIVDB protein sequences:
protease_hivdb="PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMNLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNF"

rt_hivdb = "PISPIETVPVKLKPGMDGPKVKQWPLTEEKIKALVEICTEMEKEGKISKI"\
     "GPENPYNTPVFAIKKKDSTKWRKLVDFRELNKRTQDFWEVQLGIPHPAGL"\
     "KKKKSVTVLDVGDAYFSVPLDKDFRKYTAFTIPSINNETPGIRYQYNVLP"\
     "QGWKGSPAIFQSSMTKILEPFRKQNPDIVIYQYMDDLYVGSDLEIGQHRT"\
     "KIEELRQHLLRWGFTTPDKKHQKEPPFLWMGYELHPDKWTVQPIVLPEKD"\
     "SWTVNDIQKLVGKLNWASQIYAGIKVKQLCKLLRGTKALTEVIPLTEEAE"\
     "LELAENREILKEPVHGVYYDPSKDLIAEIQKQGQGQWTYQIYQEPFKNLK"\
     "TGKYARMRGAHTNDVKQLTEAVQKIATESIVIWGKTPKFKLPIQKETWEA"\
     "WWTEYWQATWIPEWEFVNTPPLVKLWYQLEKEPIVGAETFYVDGAANRET"\
     "KLGKAGYVTDRGRQKVVSLTDTTNQKTELQAIHLALQDSGLEVNIVTDSQ"\
     "YALGIIQAQPDKSESELVSQIIEQLIKKEKVYLAWVPAHKGIGGNEQVDK"\
     "LVSAGIRKVL"

integrase_hivdb = "FLDGIDKAQEEHEKYHSNWRAMASDFNLPPVVAKEIVASCDKCQLKGEAM"\
            "HGQVDCSPGIWQLDCTHLEGKIILVAVHVASGYIEAEVIPAETGQETAYF"\
            "LLKLAGRWPVKTIHTDNGSNFTSTTVKAACWWAGIKQEFGIPYNPQSQGV"\
            "VESMNKELKKIIGQVRDQAEHLKTAVQMAVFIHNFKRKGGIGGYSAGERI"\
            "VDIIATDIQTKELQKQITKIQNFRVYYRDSRDPLWKGPAKLLWKGEGAVV"\
            "IQDNSDIKVVPRRKAKIIRDYGKQMAGDDCVASRQDED"

# #HXB2 translated protein sequences (translated based on gene coordinates above)
#protease_hxb2="PQVTLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNF"
#rt_hxb2= "PISPIETVPVKLKPGMDGPKVKQWPLTEEKIKALVEICTEMEKEGKISKI"\
#          "GPENPYNTPVFAIKKKDSTKWRKLVDFRELNKRTQDFWEVQLGIPHPAGL"\
#          "KKKKSVTVLDVGDAYFSVPLDEDFRKYTAFTIPSINNETPGIRYQYNVLP"\
#          "QGWKGSPAIFQSSMTKILEPFRKQNPDIVIYQYMDDLYVGSDLEIGQHRT"\
#          "KIEELRQHLLRWGLTTPDKKHQKEPPFLWMGYELHPDKWTVQPIVLPEKD"\
#          "SWTVNDIQKLVGKLNWASQIYPGIKVRQLCKLLRGTKALTEVIPLTEEAE"\
#          "LELAENREILKEPVHGVYYDPSKDLIAEIQKQGQGQWTYQIYQEPFKNLK"\
#          "TGKYARMRGAHTNDVKQLTEAVQKITTESIVIWGKTPKFKLPIQKETWET"\
#          "WWTEYWQATWIPEWEFVNTPPLVKLWYQLEKEPIVGAETFYVDGAANRET"\
#          "KLGKAGYVTNRGRQKVVTLTDTTNQKTELQAIYLALQDSGLEVNIVTDSQ"\
#          "YALGIIQAQPDQSESELVNQIIEQLIKKEKVYLAWVPAHKGIGGNEQVDK"\
#          "LVSAGIRKVL"
#integrase_hxb2 = "FLDGIDKAQDEHEKYHSNWRAMASDFNLPPVVAKEIVASCDKCQLKGEAM"\
#                  "HGQVDCSPGIWQLDCTHLEGKVILVAVHVASGYIEAEVIPAETGQETAYF"\
#                  "LLKLAGRWPVKTIHTDNGSNFTGATVRAACWWAGIKQEFGIPYNPQSQGV"\
#                  "VESMNKELKKIIGQVRDQAEHLKTAVQMAVFIHNFKRKGGIGGYSAGERI"\
#                  "VDIIATDIQTKELQKQITKIQNFRVYYRDSRNPLWKGPAKLLWKGEGAVV"\
#                  "IQDNSDIKVVPRRKAKIIRDYGKQMAGDDCVASRQDED"

protease_hxb2 = "CCTCAGGTCACTCTTTGGCAACGACCCCTCGTCACAATAAAGATAGGGGGGCAACTAAAGGAAGCTCTATTAGATACAGGAGCAGATGATACAGTATTAGAAGAAATGAGTTTGCCAGGAAGATGGAAACCAAAAATGATAGGGGGAATTGGAGGTTTTATCAAAGTAAGACAGTATGATCAGATACTCATAGAAATCTGTGGACATAAAGCTATAGGTACAGTATTAGTAGGACCTACACCTGTCAACATAATTGGAAGAAATCTGTTGACTCAGATTGGTTGCACTTTAAATTTT"
rt_hxb2 = "CCCATTAGCCCTATTGAGACTGTACCAGTAAAATTAAAGCCAGGAATGGATGGCCCAAAAGTTAAACAATGGCCATTGACAGAAGAAAAAATAAAAGCATTAGTAGAAATTTGTACAGAGATGGAAAAGGAAGGGAAAATTTCAAAAATTGGGCCTGAAAATCCATACAATACTCCAGTATTTGCCATAAAGAAAAAAGACAGTACTAAATGGAGAAAATTAGTAGATTTCAGAGAACTTAATAAGAGAACTCAAGACTTCTGGGAAGTTCAATTAGGAATACCACATCCCGCAGGGTTAAAAAAGAAAAAATCAGTAACAGTACTGGATGTGGGTGATGCATATTTTTCAGTTCCCTTAGATGAAGACTTCAGGAAGTATACTGCATTTACCATACCTAGTATAAACAATGAGACACCAGGGATTAGATATCAGTACAATGTGCTTCCACAGGGATGGAAAGGATCACCAGCAATATTCCAAAGTAGCATGACAAAAATCTTAGAGCCTTTTAGAAAACAAAATCCAGACATAGTTATCTATCAATACATGGATGATTTGTATGTAGGATCTGACTTAGAAATAGGGCAGCATAGAACAAAAATAGAGGAGCTGAGACAACATCTGTTGAGGTGGGGACTTACCACACCAGACAAAAAACATCAGAAAGAACCTCCATTCCTTTGGATGGGTTATGAACTCCATCCTGATAAATGGACAGTACAGCCTATAGTGCTGCCAGAAAAAGACAGCTGGACTGTCAATGACATACAGAAGTTAGTGGGGAAATTGAATTGGGCAAGTCAGATTTACCCAGGGATTAAAGTAAGGCAATTATGTAAACTCCTTAGAGGAACCAAAGCACTAACAGAAGTAATACCACTAACAGAAGAAGCAGAGCTAGAACTGGCAGAAAACAGAGAGATTCTAAAAGAACCAGTACATGGAGTGTATTATGACCCATCAAAAGACTTAATAGCAGAAATACAGAAGCAGGGGCAAGGCCAATGGACATATCAAATTTATCAAGAGCCATTTAAAAATCTGAAAACAGGAAAATATGCAAGAATGAGGGGTGCCCACACTAATGATGTAAAACAATTAACAGAGGCAGTGCAAAAAATAACCACAGAAAGCATAGTAATATGGGGAAAGACTCCTAAATTTAAACTGCCCATACAAAAGGAAACATGGGAAACATGGTGGACAGAGTATTGGCAAGCCACCTGGATTCCTGAGTGGGAGTTTGTTAATACCCCTCCCTTAGTGAAATTATGGTACCAGTTAGAGAAAGAACCCATAGTAGGAGCAGAAACCTTCTATGTAGATGGGGCAGCTAACAGGGAGACTAAATTAGGAAAAGCAGGATATGTTACTAATAGAGGAAGACAAAAAGTTGTCACCCTAACTGACACAACAAATCAGAAGACTGAGTTACAAGCAATTTATCTAGCTTTGCAGGATTCGGGATTAGAAGTAAACATAGTAACAGACTCACAATATGCATTAGGAATCATTCAAGCACAACCAGATCAAAGTGAATCAGAGTTAGTCAATCAAATAATAGAGCAGTTAATAAAAAAGGAAAAGGTCTATCTGGCATGGGTACCAGCACACAAAGGAATTGGAGGAAATGAACAAGTAGATAAATTAGTCAGTGCTGGAATCAGGAAAGTACTA"
integrase_hxb2 = "TTTTTAGATGGAATAGATAAGGCCCAAGATGAACATGAGAAATATCACAGTAATTGGAGAGCAATGGCTAGTGATTTTAACCTGCCACCTGTAGTAGCAAAAGAAATAGTAGCCAGCTGTGATAAATGTCAGCTAAAAGGAGAAGCCATGCATGGACAAGTAGACTGTAGTCCAGGAATATGGCAACTAGATTGTACACATTTAGAAGGAAAAGTTATCCTGGTAGCAGTTCATGTAGCCAGTGGATATATAGAAGCAGAAGTTATTCCAGCAGAAACAGGGCAGGAAACAGCATATTTTCTTTTAAAATTAGCAGGAAGATGGCCAGTAAAAACAATACATACTGACAATGGCAGCAATTTCACCGGTGCTACGGTTAGGGCCGCCTGTTGGTGGGCGGGAATCAAGCAGGAATTTGGAATTCCCTACAATCCCCAAAGTCAAGGAGTAGTAGAATCTATGAATAAAGAATTAAAGAAAATTATAGGACAGGTAAGAGATCAGGCTGAACATCTTAAGACAGCAGTACAAATGGCAGTATTCATCCACAATTTTAAAAGAAAAGGGGGGATTGGGGGGTACAGTGCAGGGGAAAGAATAGTAGACATAATAGCAACAGACATACAAACTAAAGAATTACAAAAACAAATTACAAAAATTCAAAATTTTCGGGTTTATTACAGGGACAGCAGAAATCCACTTTGGAAAGGACCAGCAAAGCTCCTCTGGAAAGGTGAAGGGGCAGTAGTAATACAAGATAATAGTGACATAAAAGTAGTGCCAAGAAGAAAAGCAAAGATCATTAGGGATTATGGAAAACAGATGGCAGGTGATGATTGTGTGGCAAGTAGACAGGATGAGGAT"

geneticCode ={
#Phenylalanine
"TTT": "F", "TTC": "F",
#Leucine
"TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
#Isoleucine
"ATT": "I", "ATC": "I","ATA": "I",
#Methioinine
"ATG": "M",
#Valine
"GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
#Serine
"TCT": "S", "TCC": "S", "TCA": "S","TCG": "S","AGT": "S", "AGC": "S",
#Proline
"CCT": "P", "CCC": "P", "CCA": "P","CCG": "P",
#Threonine
"ACT": "T","ACC": "T","ACA": "T","ACG": "T",
#Alanine
"GCT": "A","GCC": "A","GCA": "A","GCG": "A",
#Tyrosine
"TAT": "Y","TAC": "Y",
#Histidine
"CAT": "H","CAC": "H",
#Glutamine
"CAA": "Q","CAG": "Q",
#Asparagine
"AAT": "N","AAC": "N",
#Lysine
"AAA": "K","AAG": "K",
#Aspartic acid
"GAT": "D","GAC": "D",
#Glutamic acid
"GAA": "E","GAG": "E",
#Cysteine
"TGT": "C","TGC": "C",
#Tryptophan
"TGG": "W",
#Arginine
"CGT": "R","CGC": "R","CGA": "R","CGG": "R","AGA": "R","AGG": "R",
#Glycine
"GGT": "G","GGC": "G","GGA": "G","GGG": "G",
#stop
"TAA": "*","TAG": "*","TGA": "*"
}

def getCombinations(base1,  base2, base3,reference):
# returns all possible combinations

    mutations = []
    #print("Bases", base1, base2, base3)
    for i in range (len(base1)):
        for j in range (len(base2)):
            for k in range (len(base3)):
                codonStr = str(base1[i])+str(base2[j])+str(base3[k])
                #print("CODON:",codonStr)
                trans = geneticCode.get(codonStr,"!")
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
    ref_gene = list(refG)
    codonct = 1
    #print("LENPROT",len(prot))
    for i in range(0,len(mut),3):
        #print("I",i,prot[i])
        #print("CODONCT ",codonct-1,ref_protease[codonct-1])
        #print("Base:",prot[i],prot[i+1],prot[i+2])
        combinations = getCombinations(mut[i],mut[i+1], mut[i+2],ref_gene[codonct-1])
        if len(combinations)>0:
            # print("MUTATION:", codonct, ref_gene[codonct-1], "comb", combinations)
            mutList.append([codonct,ref_gene[codonct-1],combinations])

        codonct = codonct + 1

    #print("M",mutList)
    return mutList


def createHIVDBRequest(prot,rt,integrase, out, aavf):
    cmd = "sierrapy mutations "
    #cmd = ""
    #protease
    # print("PROT1 ",len(prot))
    
    for i in prot:
        # print("PROT: ", i)
        for j in i[2]:
            #print("SIERRA: " , i[0],i[1], j)
            mut = "PR:"+i[1]+str(i[0])+j
            #print("SIERRA", mut)
            cmd = cmd + " "+mut
            aavf.write("hxb2\t")
            aavf.write("PR\t"+str(i[0])+"\t"+i[1]+"\t"+j)
            #TODO: write filter, alt_freq, coverage, and info
            #filter: filter status; PASS if this position has passed all filters i.e. a call is made at this position. Otherwise, ...
            #alt_freq: alternate amino acid frequency; frequency of the alternate allele (float, required)
            #coverage: coverage at that position; number of reads that cover the POS (int, required)
            aavf.write("\n")
    #RT
    for i in rt:
        #print("PROT: ", i)
        for j in i[2]:
        #print("SIERRA: " , i[0],i[1], j)
            mut = "RT:"+i[1]+str(i[0])+j
            #print("SIERRA", mut)
            cmd = cmd + " "+mut
            aavf.write("hxb2\t")
            aavf.write("RT\t"+str(i[0])+"\t"+i[1]+"\t"+j)
            #TODO: write filter, alt_freq, coverage, and info
            aavf.write("\n")
    #RT
    for i in integrase:
        #print("PROT: ", i)
        for j in i[2]:
            #print("SIERRA: " , i[0],i[1], j)
            mut = "IN:"+i[1]+str(i[0])+j
            #print("SIERRA", mut)
            cmd = cmd + " "+mut
            aavf.write("hxb2\t")
            aavf.write("IN\t"+str(i[0])+"\t"+i[1]+"\t"+j)
            #TODO: write filter, alt_freq, coverage, and info
            aavf.write("\n")
    #protStr
    #cmd = cmd + " -o "+out
    #cmd = 'sierrapy mutations PR:L10I -o output_171.json'
    # print(cmd)
    #print (subprocess.check_output(cmd, shell=True))
    #output = subprocess.Popen(cmd, shell=True)

    output = subprocess.check_output(cmd, shell=True)
    file = open(out, 'w')
    file.write(output.decode("utf-8"))
    file.close()
    aavf.close()
    print(cmd)

    os.system(cmd + " -o " + out)

if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--vcfFile', required=True)
    parser.add_argument('--qual', required=False, type=int, default=None)
    parser.add_argument('--outFile', required=True)
    args = parser.parse_args()


    vcf_reader = vcf.Reader(open(args.vcfFile,'r'))

    aavffile = args.outFile.replace(".json",".aavf")
    aavf = open(aavffile, 'w')

    aavf.write("##fileformat=AAVFv1.0\n")
    aavf.write("##fileDate="+datetime.datetime.now().strftime("%Y%m%d")+"\n")
    aavf.write("##source=vcfToHivdb.py\n")
    # aavf.write("##reference="+vcf_reader.metadata['reference']+"\n")
    # for info in vcf_reader.infos:
    #         header_line = vcf_reader.infos[info].header
    #         aavf.write(header_line + '\n')
    # for filters in vcf_reader.filters:
    #         header_line = vcf_reader.filters[filters].header
    #         aavf.write(header_line + '\n')
    aavf.write("#CHROM\tGENE\tPOS\tREF\tALT\tFILTER\tALT_FREQ\tCOVERAGE\tINFO\n")

    # vcf_reader = VariantFile(snakemake.input[0])
    protease_seq = list(protease_hxb2)
    rt_seq = list(rt_hxb2)
    i_seq = list(integrase_hxb2)

    if(args.qual != None):
        for record in vcf_reader:
            #print(record)
            if record.POS >=2253 and record.POS <2550: #protease
                alt = record.ALT
                temp = []
                for i in range(len(alt[0])):
                    base = str(alt[0])[i]
                    if base not in temp:
                        temp.append(base)
                if record.QUAL <= args.qual:
                    protease_seq[record.POS-2253] = temp
            elif record.POS >=2550 and record.POS < 4229: #3869: #reverse transcriptase
                alt = record.ALT
                temp = []
                for i in range(len(alt[0])):
                    base = str(alt[0])[i]
                    if base not in temp:
                        temp.append(base)
                if record.QUAL <= args.qual:
                    rt_seq[record.POS-2550] = temp
            elif record.POS >=4230 and record.POS <=5093: #5096 minus stop codon #integrase
                alt = record.ALT
                temp = []
                for i in range(len(alt[0])):
                    base = str(alt[0])[i]
                    if base not in temp:
                        temp.append(base)
                if record.QUAL <= args.qual:
                    i_seq[record.POS-4230] = temp
    else:
        for record in vcf_reader:
            #print(record)
            if record.POS >=2253 and record.POS <2550: #protease
                alt = record.ALT
                temp = []
                for i in range(len(alt[0])):
                    base = str(alt[0])[i]
                    if base not in temp:
                        temp.append(base)
                protease_seq[record.POS-2253] = temp
            elif record.POS >=2550 and record.POS < 4229: #3869: #reverse transcriptase
                alt = record.ALT
                temp = []
                for i in range(len(alt[0])):
                    base = str(alt[0])[i]
                    if base not in temp:
                        temp.append(base)
                rt_seq[record.POS-2550] = temp
            elif record.POS >=4230 and record.POS <=5093: #5096 minus stop codon #integrase
                alt = record.ALT
                temp = []
                for i in range(len(alt[0])):
                    base = str(alt[0])[i]
                    if base not in temp:
                        temp.append(base)
                i_seq[record.POS-4230] = temp

    protease_mut = getMutations_mult(protease_seq,protease_hivdb)
    rt_mut = getMutations_mult(rt_seq,rt_hivdb)
    int_mut = getMutations_mult(i_seq,integrase_hivdb)

    #createHIVDBRequest(protease_mut, rt_mut,int_mut,snakemake.output[0])
    createHIVDBRequest(protease_mut, rt_mut, int_mut, args.outFile, aavf)




# def getHivdbPos(geneStart,vcfPos):
#     return (vcfPos-geneStart)//3
#
# def getMutAA(geneStart,vcfPos):

#
# proteaseMutations = []
# rtMutations = []
# integraseMutatoins = []

# for vcfPos: #can only look at substitutions
#     if pStart <= vcfPos <= pEnd:
#         originalAminoAcid = protease_hivdb[getHivdbPos(pStart,vcfPos)]
#         mutatedAminoAcid = getMutAA(pStart,vcfPos)
#     elif:
