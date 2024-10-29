import os
import pandas as pd
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed


#read in Goodman data
columns = ["CDS.seq", "Prot.FCC"]
df = pd.read_csv('1241934tables1.csv', usecols = columns)
df = df.rename(columns={"Prot.FCC": "Prot_FCC"})
df1 = df.dropna()
df1['CDS.seq'] = df1['CDS.seq'].str.lower()

outfile_head = ""
for k in range(2,12):
	outfile_head = f"{outfile_head},pos{k}"
	

file_object1 = open("GC_by_codon_Goodman.csv", "w")
file_object1.write(f"sequence,protein_level{outfile_head}\n")
	
for index, row in df1.iterrows():
	s= row["CDS.seq"]
	p = row["Prot_FCC"]
	#print(s,p)
	dict_gc ={}
	codons = [s[j:j+3] for j in range(3, 33, 3)]
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
