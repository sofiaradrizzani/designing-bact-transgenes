import pandas as pd
import math
import statistics

#define which columns of data i want from the source csv file
#create new dataframe using these columns
columns = ["five_prime_seq"]
df = pd.read_csv("Ecoli_five.csv", usecols = columns)

df = df.rename(columns={"five_prime_seq": "five_prime_seq"})

#extract top and bottom 25% of data points based on Protein/RNA

## dropped codon analysis
df = df[df['five_prime_seq'].notna()]
seqs = df["five_prime_seq"]

file_object2 = open("native_LO_drop.csv", "w")
file_object2.write("aa,cdn,mean_LOR_not,sd,N,sem,ends_with\n")

#create list of amino acids and dictionary of synonymous codon blocks
amino_acids = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "N", "P", "Q", "R", "S", "T", "V", "Y"]

codons_dict = {
"A": ["GCT", "GCC", "GCA", "GCG"],
"C": ["TGT", "TGC"],
"D": ["GAT", "GAC"],
"E": ["GAA", "GAG"],
"F": ["TTT", "TTC"],
"G": ["GGT", "GGC", "GGA", "GGG"],
"H": ["CAT", "CAC"],
"I": ["ATA", "ATT", "ATC"],
"K": ["AAA", "AAG"],
"L": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"],
"N": ["AAT", "AAC"],
"P": ["CCT", "CCC", "CCA", "CCG"],
"Q": ["CAA", "CAG"],
"R": ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
"S": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
"T": ["ACT", "ACC", "ACA", "ACG"],
"V": ["GTT", "GTC", "GTA", "GTG"],
"Y": ["TAT", "TAC"]
}

LOR_dropped_mn = {}
LOR_dropped_sd = {}
LOR_dropped_N = {}

for aa in amino_acids:
	codons = codons_dict[aa]
	alpha = codons
	for cdn_f in alpha:

		dropped_vals_mean_all = []
		for seq in seqs:
			seq = seq.upper()
			not_vals = []
			n = 3
			codons_high = [seq[i:i + n] for i in range(0, 30, n)]

			if cdn_f in codons_high:
				while cdn_f in codons_high:
					codons_high.remove(cdn_f)

				for cdn in codons_high:
					if cdn in dict_LOR:
						not_vals.append(dict_LOR[cdn])
				mn_not = statistics.mean(not_vals)

				dropped_vals_mean_all.append(mn_not)

		mn_codon = statistics.mean(dropped_vals_mean_all)
		sd_codon = statistics.stdev(dropped_vals_mean_all)
		N_codon = len(dropped_vals_mean_all)

		LOR_dropped_mn[cdn_f] = mn_codon
		LOR_dropped_sd[cdn_f] = sd_codon
		LOR_dropped_N[cdn_f] = N_codon
		ends_with = cdn_f[2]
		file_object2.write(
			f"{aa},{cdn_f},{mn_codon},{sd_codon},{N_codon},{sd_codon / math.sqrt(N_codon)},{ends_with}\n")

file_object2.close()

