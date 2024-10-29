import pandas as pd
import math
import os



cwd = os.getcwd()
df = pd.read_csv("Goodman_split.csv")

#extract top and bottom 25% of data points based on Protein
q75 = df["trans"].quantile(q=0.75)
high = df[df["trans"].ge(q75)]

q25 = df["trans"].quantile(q=0.25)
low = df[df["trans"].le(q25)]


#analyse each codon position
positions = ["cod2", "cod3", "cod4", "cod5", "cod6", "cod7", "cod8", "cod9", "cod10", "cod11"]

for pos in positions:
	# declare empty lists for high and low expressed codons
	codons_low_list = []
	codons_high_list = []

	# extract sequences (codons) for bottom 25% of data points, copy to codons_low_list
	# extract sequences (codons) for top 25% of data points, copy to codons_high_list
	sequences_low = low[pos].to_list()
	codons_low_list = sequences_low

	sequences_high = high[pos].to_list()
	codons_high_list = sequences_high

	# open new file to output log odds ratio for each codon
	outfile_path = os.path.join(cwd, "LO_split", "Goodman")
	os.makedirs(outfile_path, exist_ok=True)
	outfilename = f"{outfile_path}/logodds_HvL5_{pos}.csv"
	outfile = open(outfilename, "w")
	outfile.write("codon,amino_acid,log_odds,ends_with,standard_error\n")

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

	# calculate log odds ratio for the enrichment of each codon in the highly expressed group compared to synonyms
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

			if (c_high == 0) and (c_low == 0) and (not_c_high == 0) and (not_c_low == 0):
				odds_ratio = "NA"
				log_odds = "NA"
				standard_error = "NA"
			else:
				if (c_high > 0) and (c_low > 0) and (not_c_high > 0) and (not_c_low > 0):
					c_high = c_high
					c_low = c_low
					not_c_high = not_c_high
					not_c_low = not_c_low
				else:
					c_high = c_high + 0.5
					c_low = c_low + 0.5
					not_c_high = not_c_high + 0.5
					not_c_low = not_c_low + 0.5

				odds_ratio = (c_high / c_low) / (not_c_high / not_c_low)
				log_odds = round((math.log(odds_ratio)), 3)
				standard_error = math.sqrt((1 / c_high) + (1 / c_low) + (1 / not_c_high) + (1 / not_c_low))

			ends_with = c[2]

			outfile.write(f"{c},{aa},{log_odds},{ends_with},{standard_error}\n")

	outfile.close()


