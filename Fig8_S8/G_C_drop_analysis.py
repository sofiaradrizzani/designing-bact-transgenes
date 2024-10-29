import pandas as pd
import math
import statistics

#define which columns of data i want from the source csv file
#create new dataframe using these columns
columns = ["CDS.seq"]
df = pd.read_csv("1241934tables1.csv", usecols = columns)

columns1 = ["codon", "log_odds"]
df1 = pd.read_csv("Cambray_LOstab.csv", usecols = columns1)

dict_df1 =df1.set_index('codon')['log_odds'].to_dict()



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


df = df[df['CDS.seq'].notna()]
seqs = df["CDS.seq"]

file_object2 = open("Cambray_Goodman_LOstab_drop.csv", "w")
file_object2.write("aa,cdn,mean_LOR_not,sd,N,sem,ends_with\n")

LOR_dropped_mn ={}
LOR_dropped_sd ={}
LOR_dropped_N ={}

for aa in amino_acids:
	codons = codons_dict[aa]
	alpha = codons

	for cdn_f in alpha:
		dropped_vals_mean_all = []

		for seq in seqs:
			seq =seq.upper()
			not_vals= []
			n = 3
			codons_high = [seq[i:i+n] for i in range(0, 33, n)]
			
			if cdn_f in codons_high:
				codons_high.remove(cdn_f)

				for cdn in codons_high:
					if cdn in dict_df1:
						not_vals.append(dict_df1[cdn])

				mn_not = statistics.mean(not_vals)
		
				dropped_vals_mean_all.append(mn_not)
			
		mn_codon = statistics.mean(dropped_vals_mean_all)
		sd_codon = statistics.stdev(dropped_vals_mean_all)
		N_codon = len(dropped_vals_mean_all)
		
		LOR_dropped_mn[cdn_f]=mn_codon
		LOR_dropped_sd[cdn_f]=sd_codon
		LOR_dropped_N[cdn_f]=N_codon
		ends_with = cdn_f[2]	
		file_object2.write(f"{aa},{cdn_f},{mn_codon},{sd_codon},{N_codon},{sd_codon/math.sqrt(N_codon)},{ends_with}\n")
	
file_object2.close()	

		



