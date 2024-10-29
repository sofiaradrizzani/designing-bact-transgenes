import os
import pandas as pd
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed


#read in Cambray data

columns = ["gs.sequence", "clean.lin.prot.mean", "ss.rna.dna.mean"]
df = pd.read_csv("Cambray_DataS15.csv", usecols = columns)


df = df.rename(columns={"clean.lin.prot.mean": "clean_lin_prot_mean","gs.sequence": "CDS.seq", "ss.rna.dna.mean": "ss_rna_dna_mean"})

df['clean_lin_prot_mean_norm'] = df['clean_lin_prot_mean']/df['ss_rna_dna_mean']


df1 = df.dropna()
df1['CDS.seq'] = df1['CDS.seq'].str.lower()

outfile_head = ""
for k in range(2,12):
	outfile_head = f"{outfile_head},pos{k}"
	

file_object1 = open("GC_by_codon_Cambray.csv", "w")
file_object1.write(f"sequence,protein_level{outfile_head}\n")
	
for index, row in df1.iterrows():
	s= row["CDS.seq"]
	p = row['clean_lin_prot_mean_norm']
	#print(s,p)
	dict_gc ={}
	codons = [s[j:j+3] for j in range(0, 30, 3)]
	num = 2
	for cdn in codons:
		dict_gc[num] = (cdn.count('g') + cdn.count('c'))/len(cdn)		
		num = num +1
	outfile_line = f"{s},{p}"
	for k in range(2,12):
		outfile_line = f"{outfile_line},{dict_gc[k]}"
	outfile_line = f"{outfile_line}\n"
	file_object1.write(f"{outfile_line}")
file_object1.close()			
