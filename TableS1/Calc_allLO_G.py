import pandas as pd
import math

# define which columns of data i want from the source csv file
# create new dataframe using these columns
columns = ["CDS.seq", "Prot.FCC", "Trans", "dG" ,"dG.noutr", "dG.unif"]
df = pd.read_csv("1241934tables1.csv", usecols = columns)

df = df.rename(columns={"Prot.FCC": "Prot.FCC","CDS.seq": "CDS.seq", "Trans": "Trans", "dG":"dG", "dG.noutr":"dG.noutr", "dG.unif":"dG.unif"})


## dG
# extract top and bottom 25% of data points based on stability measure dG
# df['clean_lin_prot_mean_norm'] = df['clean_lin_prot_mean']/df['ss_rna_dna_mean']
q75 = df["dG"].quantile(q = 0.75)
high = df[df["dG"].ge(q75)]

q25 = df["dG"].quantile(q = 0.25)
low = df[df["dG"].le(q25)]

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
file_object1 = open("Goodman_LOstab_dG.csv", "w")
file_object1.write("codon,amino_acid,ends_with,LO_dG\n")


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

		dict_LOR[c] = log_odds

		file_object1.write(f"{c},{aa},{ends_with},{log_odds}\n")

file_object1.close()



##dG.noutr
# extract top and bottom 25% of data points based on stability measure dG.noutr
# df['clean_lin_prot_mean_norm'] = df['clean_lin_prot_mean']/df['ss_rna_dna_mean']
q75 = df["dG.noutr"].quantile(q = 0.75)
high = df[df["dG.noutr"].ge(q75)]

q25 = df["dG.noutr"].quantile(q = 0.25)
low = df[df["dG.noutr"].le(q25)]

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
file_object2 = open("Goodman_LOstab_dGnoutr.csv", "w")
file_object2.write("codon,amino_acid,ends_with,LO_dGnoutr\n")

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

		dict_LOR[c] = log_odds

		file_object2.write(f"{c},{aa},{ends_with},{log_odds}\n")

file_object2.close()



##dG.unif
# extract top and bottom 25% of data points based on stability measure dG.unif
# df['clean_lin_prot_mean_norm'] = df['clean_lin_prot_mean']/df['ss_rna_dna_mean']
q75 = df["dG.unif"].quantile(q = 0.75)
high = df[df["dG.unif"].ge(q75)]

q25 = df["dG.unif"].quantile(q = 0.25)
low = df[df["dG.unif"].le(q25)]

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
file_object3 = open("Goodman_LOstab_dGunif.csv", "w")
file_object3.write("codon,amino_acid,ends_with,LO_dGunif\n")

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

		dict_LOR[c] = log_odds

		file_object3.write(f"{c},{aa},{ends_with},{log_odds}\n")

file_object3.close()


## Prot.FCC
# extract top and bottom 25% of data points based on stability measure dG
# df['clean_lin_prot_mean_norm'] = df['clean_lin_prot_mean']/df['ss_rna_dna_mean']
q75 = df["Prot.FCC"].quantile(q = 0.75)
high = df[df["Prot.FCC"].ge(q75)]

q25 = df["Prot.FCC"].quantile(q = 0.25)
low = df[df["Prot.FCC"].le(q25)]

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
file_object4 = open("Goodman_LO_protFCC.csv", "w")
file_object4.write("codon,amino_acid,ends_with,LO_protFCC\n")


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

		dict_LOR[c] = log_odds

		file_object4.write(f"{c},{aa},{ends_with},{log_odds}\n")

file_object4.close()



## Trans
# extract top and bottom 25% of data points based on stability measure dG
# df['clean_lin_prot_mean_norm'] = df['clean_lin_prot_mean']/df['ss_rna_dna_mean']
q75 = df["Trans"].quantile(q = 0.75)
high = df[df["Trans"].ge(q75)]

q25 = df["Trans"].quantile(q = 0.25)
low = df[df["Trans"].le(q25)]

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
file_object5 = open("Goodman_LO_trans.csv", "w")
file_object5.write("codon,amino_acid,ends_with,LO_trans\n")


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

		dict_LOR[c] = log_odds

		file_object5.write(f"{c},{aa},{ends_with},{log_odds}\n")

file_object5.close()


## calculated RNAfold stability data
#read in file
df2 = pd.read_csv("GSTR.csv")

# extract top and bottom 25% of data points based on stability
# df['clean_lin_prot_mean_norm'] = df['clean_lin_prot_mean']/df['ss_rna_dna_mean']
q75 = df2["STR_5_G"].quantile(q=0.75)
high = df2[df2["STR_5_G"].ge(q75)]

q25 = df2["STR_5_G"].quantile(q=0.25)
low = df2[df2["STR_5_G"].le(q25)]

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
file_object6 = open("Goodman_LOstab_RNAfold.csv", "w")
file_object6.write("codon,amino_acid,ends_with,LO_RNAfold_G\n")


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

		dict_LOR[c] = log_odds

		file_object6.write(f"{c},{aa},{ends_with},{log_odds}\n")

file_object6.close()




