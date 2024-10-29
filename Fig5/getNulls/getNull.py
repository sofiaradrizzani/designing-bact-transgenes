import pandas as pd
from pandas import Series
import math
import string
import statistics





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



#get the Cambray log odds ratios
columns = ["codon","log_odds"]
df1 = pd.read_csv("Cambray_logodds_high_low.csv", usecols = columns)

dict_LOR =df1.set_index('codon')['log_odds'].to_dict()



file_object2 = open(f"Cambray_null.csv", "w")
file_object2.write("codon,null\n")


null_values = []	
for aa in amino_acids:
	
	codons = codons_dict[aa]
	alpha = codons
	for cdn_f in alpha:
		r = dict_LOR.copy()
		del r[cdn_f] 
		
		mean_null = Series([r[k] for k in r]).mean()
		null_values.append(mean_null)	
		
		file_object2.write(f"{cdn_f},{mean_null}\n")
	
	

file_object2.close()

dict_LOR ={}

sd_null = statistics.stdev(null_values)
mean_null = statistics.mean(null_values)
print(f'sd of Cambray null vector = {sd_null}, mean = {mean_null}')	

#get the Cambray log odds ratios
columns = ["codon","log_odds"]
df1 = pd.read_csv("Goodman_logodds_high_low.csv", usecols = columns)

dict_LOR =df1.set_index('codon')['log_odds'].to_dict()



file_object2 = open(f"Goodman_null.csv", "w")
file_object2.write("codon,null\n")


null_values = []	
for aa in amino_acids:
	
	codons = codons_dict[aa]
	alpha = codons
	for cdn_f in alpha:
		r = dict_LOR.copy()
		del r[cdn_f] 
		
		mean_null = Series([r[k] for k in r]).mean()
		null_values.append(mean_null)	
		
		file_object2.write(f"{cdn_f},{mean_null}\n")
	
	

file_object2.close()

sd_null = statistics.stdev(null_values)
mean_null = statistics.mean(null_values)
print(f'sd of Goodman null vector = {sd_null}, mean = {mean_null}')		

	



