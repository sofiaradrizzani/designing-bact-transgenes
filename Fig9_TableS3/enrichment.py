import glob
import os
import re
from copy import deepcopy
import math

cwd = os.getcwd()

# labelling amino acids (aas) to codons and making a table of all codons
cdns_I = ["ATA", "ATC", "ATT"]
cdns_T = ["ACA", "ACC", "ACG", "ACT"]
cdns_N = ['AAC', 'AAT']
cdns_K = ['AAA', 'AAG']
cdns_S = ['AGC', 'AGT', 'TCA', 'TCC', 'TCG', 'TCT']
cdns_R = ['AGA', 'AGG', 'CGA', 'CGC', 'CGG', 'CGT']
cdns_L = ['CTA', 'CTC', 'CTG', 'CTT', 'TTA', 'TTG']
cdns_P = ['CCA', 'CCC', 'CCG', 'CCT']
cdns_H = ['CAC', 'CAT']
cdns_Q = ['CAA', 'CAG']
cdns_V = ['GTA', 'GTC', 'GTG', 'GTT']
cdns_A = ['GCA', 'GCC', 'GCG', 'GCT']
cdns_D = ['GAC', 'GAT']
cdns_E = ['GAA', 'GAG']
cdns_G = ['GGA', 'GGC', 'GGG', 'GGT']
cdns_F = ['TTC', 'TTT']
cdns_Y = ['TAC', 'TAT']
cdns_C = ['TGC', 'TGT']

all_cdns = [cdns_I, cdns_T, cdns_N, cdns_K, cdns_S, cdns_R, cdns_L, cdns_P, cdns_H, cdns_Q, cdns_V, cdns_A, cdns_D,
            cdns_E, cdns_G, cdns_F, cdns_Y, cdns_C]
aas = []

for c in all_cdns:
    aas.append(c)

table = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
    'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W',
}


os.chdir(f"{cwd}/all_genomes")

allgenera = glob.glob("*.fna")

for g in allgenera:
    genus = re.search("(.*?).fna", g)
    genus = genus.group()
    genus = re.sub(".fna", "", genus)
    print(f"Analysing {genus}")

    os.chdir(cwd)

    if not os.path.exists(f"{cwd}/enrichment"):
        os.makedirs(f"{cwd}/enrichment")
    os.chdir(f"{cwd}/enrichment")

    outfile_lo = open(f"{genus}_V5core.csv", "w")
    outfile_lo.write(f"codon,amino_acid,{genus},third_nucl,standard_error\n")

    infile_five_path = os.path.join(cwd, "all_seqs",f"{genus}_five.csv")
    infile_five = open(infile_five_path, "r")
    seqs_five = infile_five.read()
    infile_five.close()

    infile_core_path = os.path.join(cwd, "all_seqs", f"{genus}_core.csv")
    infile_core = open(infile_core_path, "r")
    seqs_core = infile_core.read()
    infile_core.close()

    parts_f = seqs_five.splitlines()
    parts_f = parts_f[1:]

    parts_c = seqs_core.splitlines()
    parts_c = parts_c[1:]


    # double-check number of genes
    #print(f"Number five primes: {len(parts_f)}")
    #print(f"Number genetic cores: {len(parts_c)}")


    count_codons_F = {}
    count_codons_C = {}
    log_odds = {}

    # get codon counts at 5'
    for g in parts_f:
        g = g.split(",")
        #seq = list(filter(lambda i: len(i) == 33, g))
        seq = g[3]
        seq = "".join(seq)
        seq = seq.upper()

        if re.search("[^GCAT].*",seq):
            print(f"{gene} rubbish sequence")
            continue
        # print(g)
        # print(seq)
        for i in range(0, len(seq), 3):
            codon = seq[i:i + 3]
            codon = codon.upper()
            if codon in count_codons_F:
                count_codons_F[codon] = count_codons_F[codon] + 1
            else:
                count_codons_F[codon] = 1
    #print(f'Five prime numbers: {count_codons_F}')

    # get codon counts at genetic core
    for g in parts_c:
        g = g.split(",")

        gene = g[0]
        prot_id = g[1]
        #seq = list(filter(lambda i: len(i) == 33, g))
        seq = g[3]
        seq = "".join(seq)
        seq = seq.upper()

        if re.search("[^GCAT].*",seq):
            print(f"{gene} rubbish sequence")
            continue
        # print(g)
        # print(seq)
        for i in range(0, len(seq), 3):
            codon = seq[i:i + 3]
            codon = codon.upper()
            if codon in count_codons_C:
                count_codons_C[codon] = count_codons_C[codon] + 1
            else:
                count_codons_C[codon] = 1
    #print(f'Core numbers: {count_codons_C}')

    # finding odds ratio and standard error
    for a in aas:
        #print(a)
        b = deepcopy(a)

        for j in b:
            if j not in count_codons_F: # some codons are not present in any sequence so need to attribute count of 0
                count_codons_F[f"{j}"] = 0

            if j not in count_codons_C:
                count_codons_C[f"{j}"] = 0

        for j in b:
            a = deepcopy(b)  # will need at a.remove step
            numj_f = count_codons_F[j]
            numj_c = count_codons_C[j]
            a.remove(j)
            #print(a)
            sum_notj_f = 0
            sum_notj_c = 0

            for cdn in a: # summing up the frequencies of the non-focal codons (not j)
                sum_notj_f = sum_notj_f + count_codons_F[cdn]
                sum_notj_c = sum_notj_c + count_codons_C[cdn]

            if (numj_f == 0) and (numj_c == 0) and (sum_notj_f == 0) and (sum_notj_c == 0):
                odds_ratio = "NA"
                log_odds = "NA"
                standard_error = "NA"
            else:
                if (numj_f > 0) and (numj_c > 0) and (sum_notj_f > 0) and (sum_notj_c > 0):
                    numj_f = numj_f
                    numj_c = numj_c
                    sum_notj_f = sum_notj_f
                    sum_notj_c = sum_notj_c
                else:
                    numj_f = numj_f + 0.5
                    numj_c = numj_c + 0.5
                    sum_notj_f = sum_notj_f + 0.5
                    sum_notj_c = sum_notj_c + 0.5

                odds_ratio = (numj_f / numj_c) / (sum_notj_f / sum_notj_c)
                log_odds = round((math.log(odds_ratio)), 3)
                standard_error = math.sqrt((1 / numj_f) + (1 / numj_c) + (1 / sum_notj_f) + (1 / sum_notj_c))

            #if (numj_f <= 0) or (numj_c <= 0):
                #odds_ratio = "NA"
                #log_odds = "NA"
                #standard_error = "NA"
            #elif (sum_notj_f <= 0) or (sum_notj_c <= 0):
                #odds_ratio = "NA"
                #log_odds = "NA"
                #standard_error = "NA"
            #else:
                #odds_ratio = (numj_f / numj_c) / (sum_notj_f / sum_notj_c)
                #log_odds = round((math.log(odds_ratio)), 3)
                #standard_error = math.sqrt((1 / numj_f) + (1 / numj_c) + (1 / sum_notj_f) + (1 / sum_notj_c))

            #print(f'{j}:{log_odds}, std error: {standard_error}')

                amino = table[j]

                third_nucl = j[2]

                outfile_lo.write(f"{j},{amino},{log_odds},{third_nucl},{standard_error}\n")

    outfile_lo.close()

    os.chdir(f"{cwd}/all_genomes")



