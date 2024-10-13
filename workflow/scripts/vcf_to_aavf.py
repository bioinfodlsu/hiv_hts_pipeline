import vcf
import argparse
import datetime

# Gene coordinates (1-based)
GENE_COORDS = {
    "PR": (2253, 2549),
    "RT": (2550, 4229),
    "IN": (4230, 5093)
}

#HIVDB protein sequences:
PROTEIN_SEQUENCES = {
    "PR": "PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMNLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNF",
    "RT": "PISPIETVPVKLKPGMDGPKVKQWPLTEEKIKALVEICTEMEKEGKISKI"\
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
            "LVSAGIRKVL",
    "IN": "FLDGIDKAQEEHEKYHSNWRAMASDFNLPPVVAKEIVASCDKCQLKGEAM"\
            "HGQVDCSPGIWQLDCTHLEGKIILVAVHVASGYIEAEVIPAETGQETAYF"\
            "LLKLAGRWPVKTIHTDNGSNFTSTTVKAACWWAGIKQEFGIPYNPQSQGV"\
            "VESMNKELKKIIGQVRDQAEHLKTAVQMAVFIHNFKRKGGIGGYSAGERI"\
            "VDIIATDIQTKELQKQITKIQNFRVYYRDSRDPLWKGPAKLLWKGEGAVV"\
            "IQDNSDIKVVPRRKAKIIRDYGKQMAGDDCVASRQDED"
}

# HIVDB reference sequences
REFERECE_SEQUENCES = {
    "PR": "CCTCAGGTCACTCTTTGGCAACGACCCCTCGTCACAATAAAGATAGGGGGGCAACTAAAGGAAGCTCTATTAGATACAGGAGCAGATGATACAGTATTAGAAGAAATGAGTTTGCCAGGAAGATGGAAACCAAAAATGATAGGGGGAATTGGAGGTTTTATCAAAGTAAGACAGTATGATCAGATACTCATAGAAATCTGTGGACATAAAGCTATAGGTACAGTATTAGTAGGACCTACACCTGTCAACATAATTGGAAGAAATCTGTTGACTCAGATTGGTTGCACTTTAAATTTT",
    "RT": "CCCATTAGCCCTATTGAGACTGTACCAGTAAAATTAAAGCCAGGAATGGATGGCCCAAAAGTTAAACAATGGCCATTGACAGAAGAAAAAATAAAAGCATTAGTAGAAATTTGTACAGAGATGGAAAAGGAAGGGAAAATTTCAAAAATTGGGCCTGAAAATCCATACAATACTCCAGTATTTGCCATAAAGAAAAAAGACAGTACTAAATGGAGAAAATTAGTAGATTTCAGAGAACTTAATAAGAGAACTCAAGACTTCTGGGAAGTTCAATTAGGAATACCACATCCCGCAGGGTTAAAAAAGAAAAAATCAGTAACAGTACTGGATGTGGGTGATGCATATTTTTCAGTTCCCTTAGATGAAGACTTCAGGAAGTATACTGCATTTACCATACCTAGTATAAACAATGAGACACCAGGGATTAGATATCAGTACAATGTGCTTCCACAGGGATGGAAAGGATCACCAGCAATATTCCAAAGTAGCATGACAAAAATCTTAGAGCCTTTTAGAAAACAAAATCCAGACATAGTTATCTATCAATACATGGATGATTTGTATGTAGGATCTGACTTAGAAATAGGGCAGCATAGAACAAAAATAGAGGAGCTGAGACAACATCTGTTGAGGTGGGGACTTACCACACCAGACAAAAAACATCAGAAAGAACCTCCATTCCTTTGGATGGGTTATGAACTCCATCCTGATAAATGGACAGTACAGCCTATAGTGCTGCCAGAAAAAGACAGCTGGACTGTCAATGACATACAGAAGTTAGTGGGGAAATTGAATTGGGCAAGTCAGATTTACCCAGGGATTAAAGTAAGGCAATTATGTAAACTCCTTAGAGGAACCAAAGCACTAACAGAAGTAATACCACTAACAGAAGAAGCAGAGCTAGAACTGGCAGAAAACAGAGAGATTCTAAAAGAACCAGTACATGGAGTGTATTATGACCCATCAAAAGACTTAATAGCAGAAATACAGAAGCAGGGGCAAGGCCAATGGACATATCAAATTTATCAAGAGCCATTTAAAAATCTGAAAACAGGAAAATATGCAAGAATGAGGGGTGCCCACACTAATGATGTAAAACAATTAACAGAGGCAGTGCAAAAAATAACCACAGAAAGCATAGTAATATGGGGAAAGACTCCTAAATTTAAACTGCCCATACAAAAGGAAACATGGGAAACATGGTGGACAGAGTATTGGCAAGCCACCTGGATTCCTGAGTGGGAGTTTGTTAATACCCCTCCCTTAGTGAAATTATGGTACCAGTTAGAGAAAGAACCCATAGTAGGAGCAGAAACCTTCTATGTAGATGGGGCAGCTAACAGGGAGACTAAATTAGGAAAAGCAGGATATGTTACTAATAGAGGAAGACAAAAAGTTGTCACCCTAACTGACACAACAAATCAGAAGACTGAGTTACAAGCAATTTATCTAGCTTTGCAGGATTCGGGATTAGAAGTAAACATAGTAACAGACTCACAATATGCATTAGGAATCATTCAAGCACAACCAGATCAAAGTGAATCAGAGTTAGTCAATCAAATAATAGAGCAGTTAATAAAAAAGGAAAAGGTCTATCTGGCATGGGTACCAGCACACAAAGGAATTGGAGGAAATGAACAAGTAGATAAATTAGTCAGTGCTGGAATCAGGAAAGTACTA",
    "IN": "TTTTTAGATGGAATAGATAAGGCCCAAGATGAACATGAGAAATATCACAGTAATTGGAGAGCAATGGCTAGTGATTTTAACCTGCCACCTGTAGTAGCAAAAGAAATAGTAGCCAGCTGTGATAAATGTCAGCTAAAAGGAGAAGCCATGCATGGACAAGTAGACTGTAGTCCAGGAATATGGCAACTAGATTGTACACATTTAGAAGGAAAAGTTATCCTGGTAGCAGTTCATGTAGCCAGTGGATATATAGAAGCAGAAGTTATTCCAGCAGAAACAGGGCAGGAAACAGCATATTTTCTTTTAAAATTAGCAGGAAGATGGCCAGTAAAAACAATACATACTGACAATGGCAGCAATTTCACCGGTGCTACGGTTAGGGCCGCCTGTTGGTGGGCGGGAATCAAGCAGGAATTTGGAATTCCCTACAATCCCCAAAGTCAAGGAGTAGTAGAATCTATGAATAAAGAATTAAAGAAAATTATAGGACAGGTAAGAGATCAGGCTGAACATCTTAAGACAGCAGTACAAATGGCAGTATTCATCCACAATTTTAAAAGAAAAGGGGGGATTGGGGGGTACAGTGCAGGGGAAAGAATAGTAGACATAATAGCAACAGACATACAAACTAAAGAATTACAAAAACAAATTACAAAAATTCAAAATTTTCGGGTTTATTACAGGGACAGCAGAAATCCACTTTGGAAAGGACCAGCAAAGCTCCTCTGGAAAGGTGAAGGGGCAGTAGTAATACAAGATAATAGTGACATAAAAGTAGTGCCAAGAAGAAAAGCAAAGATCATTAGGGATTATGGAAAACAGATGGCAGGTGATGATTGTGTGGCAAGTAGACAGGATGAGGAT"
}

GENETIC_CODE  ={
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

def get_combinations(base1,base2, base3,reference):
    mutations = []
    #print("Bases", base1, base2, base3)
    for i in range (len(base1)):
        for j in range (len(base2)):
            for k in range (len(base3)):
                codonStr = str(base1[i])+str(base2[j])+str(base3[k])
                trans = GENETIC_CODE.get(codonStr,"!")
                # print("CODON:",codonStr,"TRANS:",trans, reference)

                if trans!="!" and trans!=reference:
                    #add as a mutation only if it differs from the reference
                    if [trans, codonStr] not in mutations:
                        #print("ALT codon",codonStr)
                        mutations.append([trans, codonStr])
    return mutations

def parse_vcf(vcf_reader):
    gene_data = {gene: list(REFERECE_SEQUENCES[gene]) for gene in GENE_COORDS.keys()}
    record_infos = {}
    for record in vcf_reader:
        for gene, (start, end) in GENE_COORDS.items():
            if start <= record.POS <= end:
                pos_in_gene = record.POS - start
                alt = [str(alt_base) for alt_base in record.ALT if alt_base != gene_data[gene][pos_in_gene]]
                gene_data[gene][pos_in_gene] = alt
                record_infos[(gene, pos_in_gene)] = {'CHROM': record.CHROM} | record.INFO
    return gene_data, record_infos

def write_aavf_header(aavf, vcf_metadata):
    aavf.write("##fileformat=AAVFv1.0\n")
    aavf.write(f"##fileDate={datetime.datetime.now().strftime('%Y%m%d')}\n")
    aavf.write("##source=vcf_to_aavf.py\n")
    aavf.write(f"##reference={vcf_metadata['reference']}\n")
    aavf.write("##INFO=<ID=RC,Number=1,Type=String,Description=\"Reference codon\">\n")
    aavf.write("##INFO=<ID=AC,Number=.,Type=String,Description=\"Alternate codon\">\n")
    aavf.write("##INFO=<ID=ACC,Number=.,Type=Integer,Description=\"Alternate codon coverage\">\n")
    aavf.write("##INFO=<ID=ACF,Number=.,Type=Float,Description=\"Alternate codon frequency\">\n")

    # write filters here

    aavf.write("#CHROM\tGENE\tPOS\tREF\tALT\tFILTER\tALT_FREQ\tCOVERAGE\tINFO\n")

def write_mutations(aavf, gene, gene_data, record_infos, ref_prot_seq):
    codon_ct = 1
    ref_seq = list(REFERECE_SEQUENCES[gene])
    for i in range(0, len(ref_prot_seq), 3):
        ref_codon = ''.join(ref_seq[i:i+3])
        combinations = get_combinations(gene_data[i], gene_data[i+1], gene_data[i+2], ref_prot_seq[codon_ct - 1])
        for trans, alt_codon in combinations:
            record_info = [{}, {}, {}]
            for j in range(3):
                if (gene, i + j) in record_infos:
                    record_info[j] = record_infos[(gene, i + j)]
            filtered_data = [item for item in record_info if item and 'AF' in item]
            if filtered_data:
                min_af_dict = min(filtered_data, key=lambda x: x['AF'])
                chrom, alt_freq, coverage = min_af_dict['CHROM'], min_af_dict['AF'], min_af_dict['DP']
                info = f"RC={ref_codon};AC={alt_codon};ACC={coverage};ACF={alt_freq}"
                aavf.write(f"{chrom}\t{gene}\t{codon_ct}\t{ref_prot_seq[codon_ct-1]}\t{trans}\t.\t{alt_freq}\t{coverage}\t{info}\n")
        codon_ct += 1

def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--vcf', required=True, help="VCF file to convert to AAVF")
    parser.add_argument("--aavf", required=True, help="Output AAVF file")
    args = parser.parse_args()

    with open(args.vcf, 'r') as vcf_file, open(args.aavf, 'w') as aavf:
        vcf_reader = vcf.Reader(vcf_file)
        write_aavf_header(aavf, vcf_reader.metadata)
        gene_data, record_infos = parse_vcf(vcf_reader)
        for gene, ref_prot_seq in PROTEIN_SEQUENCES.items():
            write_mutations(aavf, gene, gene_data[gene], record_infos, ref_prot_seq)

if __name__ == '__main__':
    main()