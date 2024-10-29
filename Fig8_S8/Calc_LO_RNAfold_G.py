import pandas as pd
import math
import statistics


## calculated RNAfold stability data
#read in file
df = pd.read_csv("GSTR.csv")

# extract top and bottom 25% of data points based on stability
# df['clean_lin_prot_mean_norm'] = df['clean_lin_prot_mean']/df['ss_rna_dna_mean']
q75 = df["STR_5_G"].quantile(q=0.75)
high = df[df["STR_5_G"].ge(q75)]

q25 = df["STR_5_G"].quantile(q=0.25)
low = df[df["STR_5_G"].le(q25)]

# declare empty lists for high and low expressed codons
codons_low_list = []
codons_high_list = []

# extract sequences for bottom 25% of data points, turn into list, convert into codons and add to codons_low_list
# extract sequences for top 25% of data points, turn into list, convert into codons and add to codons_high_list
sequences_low = low["CDS.seq"].to_list()
for seq in sequences_low:
	seq = seq.upper()
	#print(seq)
	n = 3
	codons_low = [seq[i:i + n] for i in range(3, 33, n)]
	codons_low_list = codons_low_list + codons_low

sequences_high = high["CDS.seq"].to_list()
for seq in sequences_high:
	seq = seq.upper()
	#print(seq)
	n = 3
	codons_high = [seq[i:i + n] for i in range(3, 33, n)]
	codons_high_list = codons_high_list + codons_high

# open new file to output log odds ratio for each codon
file_object1 = open("Goodman_LOstab.csv", "w")
file_object1.write("codon,amino_acid,ends_with,log_odds\n")

# create list of amino acids and dictionary of synonymous codon blocks
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

dict_LOR = {}

# calculate log odds ratio for the enrichment of each codon in the highly stable group compared to synonyms
for aa in amino_acids:
	codons = codons_dict[aa]
	alpha = codons

	for c in alpha:
		c_high = codons_high_list.count(c)
		c_low = codons_low_list.count(c)
		not_c = [codon for codon in codons if not codon.startswith(c)]
		not_c_high = 0
		not_c_low = 0

		for nc in not_c:
			not_c_high = not_c_high + codons_high_list.count(nc)
			not_c_low = not_c_low + codons_low_list.count(nc)

		odds_ratio = (c_high / c_low) / (not_c_high / not_c_low)
		log_odds = round((math.log(odds_ratio)), 4)

		ends_with = c[2]

		dict_LOR[c] = log_odds

		file_object1.write(f"{c},{aa},{ends_with},{log_odds}\n")

file_object1.close()


## dropped codon analysis
df = df[df['CDS.seq'].notna()]
seqs = df["CDS.seq"]

file_object2 = open("Goodman_LOstab_drop.csv", "w")
file_object2.write("aa,cdn,mean_LOR_not,sd,N,sem,ends_with\n")

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

