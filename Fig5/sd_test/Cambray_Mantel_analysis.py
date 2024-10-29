import pandas as pd
import math
import string
import statistics
import random




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




#define which columns of data i want from the source csv file
#create new dataframe using these columns

#get the Cambray log odds ratios
columns = ["codon","log_odds"]
df1 = pd.read_csv("Cambray_logodds_high_low.csv", usecols = columns)

dict_LOR =df1.set_index('codon')['log_odds'].to_dict()

#get the Cambray sequences
columns = ["gs.sequence"]
df = pd.read_csv("Cambray_DataS15.csv", usecols = columns)

df = df.rename(columns={"gs.sequence": "CDS.seq"})

df = df[df['CDS.seq'].notna()]

seqs = df["CDS.seq"].tolist()

#get Goodman seq just to know sample size for drop analysis
columns = ["CDS.seq"]
df2 = pd.read_csv("1241934tables1.csv", usecols = columns)

df2 = df2.rename(columns={"CDS.seq": "CDS_seq"})

df2 = df2[df2['CDS_seq'].notna()]


Gseqs = df2["CDS_seq"]
num_seqs = len(Gseqs)
del Gseqs
del df2

#get Goodman_dropped values to calculate sd
column = ["mean_LOR_not"]
df3 = pd.read_csv("Cambray_Goodman_codondrop_LOR.csv", usecols = column)
df3 = df3[df3['mean_LOR_not'].notna()]

O_sd = statistics.stdev(df3['mean_LOR_not'])

R_sds = []

file_object2 = open(f"Cambray_codondrop_LOR_Mantel_sds.csv", "w")
file_object2.write("run,sd\n")
file_object2.write(f"Observed,{O_sd}\n")

for i in range(0,1000):
	if i % 10 ==0:
		print(f"Done {i} simulations")
	
	seqs_rand = random.sample(seqs,num_seqs)
	LOR_dropped_mn =[]
	
	for aa in amino_acids:

		codons = codons_dict[aa]
		alpha = codons
		for cdn_f in alpha:
	
			dropped_vals_mean_all = []
			for seq in seqs_rand:
				seq =seq.upper()
				not_vals= []
				n = 3
				codons_high = [seq[i:i+n] for i in range(0, 30, n)]
			
				if cdn_f in codons_high:
					while cdn_f in codons_high:
						codons_high.remove(cdn_f)


					for cdn in codons_high:
						if cdn in dict_LOR:
							not_vals.append(dict_LOR[cdn])
					mn_not = statistics.mean(not_vals)
		
					dropped_vals_mean_all.append(mn_not)
			
			mn_codon = statistics.mean(dropped_vals_mean_all)
			N_codon = len(dropped_vals_mean_all)
		
			LOR_dropped_mn.append(mn_codon)
			
			
				

			#file_object2.write(f"{aa},{cdn_f},{mn_codon},{N_codon}\n")
		
	R_sd = statistics.stdev(LOR_dropped_mn)
	R_sds.append(R_sd)
	file_object2.write(f"{i},{R_sd}\n")
	
	
	
num_simulations = len(R_sds)	
R_greaterthanobserved = [x for x in R_sds if x >= O_sd]		

P_value = len(R_greaterthanobserved)/num_simulations

Zapprox = (O_sd - statistics.mean(R_sds))/statistics.stdev(R_sds)

print(f"Observed sd = {O_sd}\nNumber of equal sized sims >= : {len(R_greaterthanobserved)} of {num_simulations}\nthus P= {P_value}, Zapprox = {Zapprox}")


file_object2.close()		



