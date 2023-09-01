import argparse
import os
import vcf
import subprocess

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


def createHIVDBRequest(prot,rt,integrase, out):
    cmd = "sierrapy mutations "
    #cmd = ""
    #protease
    for i in prot:
        # print("PROT: ", i)
        for j in i[2]:
            #print("SIERRA: " , i[0],i[1], j)
            mut = "PR:"+i[1]+str(i[0])+j
            #print("SIERRA", mut)
            cmd = cmd + " "+mut
    #RT
    for i in rt:
        #print("PROT: ", i)
        for j in i[2]:
        #print("SIERRA: " , i[0],i[1], j)
            mut = "RT:"+i[1]+str(i[0])+j
            #print("SIERRA", mut)
            cmd = cmd + " "+mut
    #I
    for i in integrase:
        #print("PROT: ", i)
        for j in i[2]:
            #print("SIERRA: " , i[0],i[1], j)
            mut = "RT:"+i[1]+str(i[0])+j
            #print("SIERRA", mut)
            cmd = cmd + " "+mut

    output = subprocess.check_output(cmd, shell=True)
    file = open(out, 'w')
    file.write(output.decode("utf-8"))
    file.close()


if __name__ == '__main__':

    # parser = argparse.ArgumentParser(description='')
    # parser.add_argument('--vcfFile', required=True)
    # parser.add_argument('--mafFile', required=True)
    # parser.add_argument('--outFile', required=True)
    # args = parser.parse_args()

    # vcf_reader = vcf.Reader(open(args.vcfFile,'r'))
    vcf_reader = vcf.Reader(open("/Users/kimwilliamelee/Git/hiv/test_out/variants/Brazil/SRR7993842/paramgroup_1/variants_last.vcf",'r'))

    # maf = open(args.mafFile, 'r')  
    file = "/Users/kimwilliamelee/Git/hiv/test_out/phase_two_subtype_alignments/hxb2_subtype02_AG_2.maf"
    maf = open(file, "r")
    mapping = []

    for line in (maf.readlines() [-3:]):
        # print(line, end ='')
        mapping.append(line.split())

    # print(mapping)


    protease_seq = list(protease_hxb2)
    rt_seq = list(rt_hxb2)
    i_seq = list(integrase_hxb2)

    for record in vcf_reader:
        #print(record)
        if record.POS >=pStart and record.POS <pEnd: #protease
            alt = record.ALT
            temp = []
            for i in range(len(alt[0])):
                base = str(alt[0])[i]
                if base not in temp:
                    temp.append(base)
            protease_seq[record.POS-2253] = temp
        elif record.POS >=rtStart and record.POS < rtEnd: #3869: #reverse transcriptase
            alt = record.ALT
            temp = []
            for i in range(len(alt[0])):
                base = str(alt[0])[i]
                if base not in temp:
                    temp.append(base)
            rt_seq[record.POS-2550] = temp
        elif record.POS >=iStart and record.POS <=iEnd: #5096 minus stop codon #integrase
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

    # createHIVDBRequest(protease_mut, rt_mut,int_mut,args.outFile)
    createHIVDBRequest(protease_mut, rt_mut,int_mut,"/Users/kimwilliamelee/Git/hiv/test_out/variants/Brazil/SRR7993842/paramgroup_1/output.txt")