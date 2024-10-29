import pandas as pd
import math
import string
import statistics
import random


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
"M": ["ATG"],
"N": ["AAT", "AAC"],
"P": ["CCT", "CCC", "CCA", "CCG"],
"Q": ["CAA", "CAG"],
"R": ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
"S": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
"T": ["ACT", "ACC", "ACA", "ACG"],
"V": ["GTT", "GTC", "GTA", "GTG"],
"W": ["TGG"],
"Y": ["TAT", "TAC"]
}

#define which columns of data i want from the source csv file
#create new dataframe using these columns

columnsG = ["CDS.seq"]
columnsC = ['gs.sequence']
dfC = pd.read_csv("Cambray_DataS15.csv", usecols = columnsC)
dfG = pd.read_csv("1241934tables1.csv", usecols = columnsG)



dfC = dfC.rename(columns={"gs.sequence": "CDS.seq"})
dfC = dfC[dfC['CDS.seq'].notna()]
dfG = dfG[dfG['CDS.seq'].notna()]



sequences_C = dfC["CDS.seq"].to_list()
sequences_G = dfG["CDS.seq"].to_list()

num_seqs = len(sequences_G)

#declare empty lists for high and low expressed codons


codons_used_G = []

count_cnds_G ={}

#see which codons are employed

			
for seq in sequences_G:
	if len(codons_used_G) ==61:
		continue
	seq =seq.upper()
	n = 3
	codons_con = [seq[i:i+n] for i in range(3, 33, n)]
	for c in codons_con:
		if c not in codons_used_G:
			codons_used_G.append(c)
for seq in sequences_G:

	seq =seq.upper()
	n = 3
	codons_con = [seq[i:i+n] for i in range(3, 33, n)]
	for c in codons_con:
		if c in count_cnds_G:
			count_cnds_G[c] += 1
		else:
			count_cnds_G[c] = 1			

total = 0
for c in codons_used_G:
	total = total + count_cnds_G[c]	

freq_codonG ={}
for c in codons_used_G:
	freq_codonG[c] = count_cnds_G[c]/total				
print(f'{len(codons_used_G)} codons employed in Goodman data')			


#make a dictionary in whcih each focal codon has a count of nonfocal

codon_cnt_droppedG ={}

constructs_with_cdnG ={}
Ecodon_cnt_droppedG={}


for cdn_f in codons_used_G:
	r = codons_used_G.copy()
	r.remove(cdn_f)
	for cdn_nonf in r:
		codon_cnt_droppedG[cdn_f,cdn_nonf] = 0
	

for cdn_f in codons_used_G:
	
	for seq in sequences_G:
		seq =seq.upper()
		not_vals= []
		n = 3
		codons_high = [seq[i:i+n] for i in range(3, 33, n)]
			
		if cdn_f in codons_high:
			if cdn_f in constructs_with_cdnG:
				constructs_with_cdnG[cdn_f] +=1
			else:
				constructs_with_cdnG[cdn_f] =1	
			
			while cdn_f in codons_high:
				codons_high.remove(cdn_f)

			for cdn in codons_high:
				codon_cnt_droppedG[cdn_f,cdn] +=1
				

	
#file_object2.close()	

#generate expected from number of codons used and number of constructs with a given codon		
#open new file to output log odds ratio for each codon



file_objectG = open("Goodman_chi.csv", "w")
file_objectG.write("codon_focal,codon_flanking,Ob,Exp,chi\n")


file_object_rand = open("Rand_chi_Efromfreq.csv", "w")
file_object_rand.write("run,chi\n")

#sum_chiG =[]
chi_valsG = []

for cdn_f in codons_used_G:
	r = codons_used_G.copy()
	r.remove(cdn_f)
	for cdn_nonf in r:
		Ecodon_cnt_droppedG[cdn_f,cdn_nonf] = (constructs_with_cdnG[cdn_f]*10 - count_cnds_G[cdn_f])*(freq_codonG[cdn_nonf]/(1-freq_codonG[cdn_f]))
		chi = (codon_cnt_droppedG[cdn_f,cdn_nonf] - Ecodon_cnt_droppedG[cdn_f,cdn_nonf])*(codon_cnt_droppedG[cdn_f,cdn_nonf] - Ecodon_cnt_droppedG[cdn_f,cdn_nonf])/Ecodon_cnt_droppedG[cdn_f,cdn_nonf]
		file_objectG.write(f"{cdn_f},{cdn_nonf},{codon_cnt_droppedG[cdn_f,cdn_nonf]},{Ecodon_cnt_droppedG[cdn_f,cdn_nonf]},{chi}\n")
		chi_valsG.append(chi)

chi_sumG =sum(chi_valsG)
file_objectG.write(f",total=,{chi_sumG}")

		
file_objectG.close()

file_object_rand.write(f"G_chi,{chi_sumG}\n")
print("Finished Goodman data")



#now subsample from Cambray data

for i in range(0,1000):
	if i % 100 ==0:
		print(f"Done {i} simulations")
	
	seqs_rand = random.sample(sequences_C,num_seqs)
	
	codons_used_C = []	
	count_cnds_C ={}
	

	for seq in seqs_rand:
		if len(codons_used_C) ==61:
			continue
		seq =seq.upper()
		n = 3
		codons_con = [seq[i:i+n] for i in range(0, 30, n)]
		for c in codons_con:
			if c not in codons_used_C:
				codons_used_C.append(c)
	for seq in seqs_rand:
		seq =seq.upper()
		n = 3
		codons_con = [seq[i:i+n] for i in range(0, 30, n)]
		for c in codons_con:
			if c in count_cnds_C:
				count_cnds_C[c] += 1
			else:
				count_cnds_C[c] = 1	
	total = 0
	for c in codons_used_C:
		total = total + count_cnds_C[c]	

	freq_codonC ={}
	for c in codons_used_C:
		freq_codonC[c] = count_cnds_C[c]/total					
	print(f'{len(codons_used_C)} codons employed in Cambray data, sim {i}')

	codon_cnt_droppedC ={}

	constructs_with_cdnC ={}

	Ecodon_cnt_droppedC={}	

	for cdn_f in codons_used_C:
		r = codons_used_C.copy()
		r.remove(cdn_f)
		for cdn_nonf in r:
			codon_cnt_droppedC[cdn_f,cdn_nonf] = 0

	for cdn_f in codons_used_C:
	
		for seq in seqs_rand:
			seq =seq.upper()
			not_vals= []
			n = 3
			codons_high = [seq[i:i+n] for i in range(0, 30, n)]
			
			if cdn_f in codons_high:
				if cdn_f in constructs_with_cdnC:
					constructs_with_cdnC[cdn_f] +=1
				else:
					constructs_with_cdnC[cdn_f] =1	
			
				while cdn_f in codons_high:
					codons_high.remove(cdn_f)

				for cdn in codons_high:
					codon_cnt_droppedC[cdn_f,cdn] +=1				

	#file_objectC = open("Cambray_chi.csv", "w")
	#file_objectC.write("codon_focal,codon_flanking,Ob,Exp,chi\n")	
	chi_valsC = []	
	for cdn_f in codons_used_C:
		r = codons_used_C.copy()
		r.remove(cdn_f)
		for cdn_nonf in r:
			Ecodon_cnt_droppedC[cdn_f,cdn_nonf] = (constructs_with_cdnC[cdn_f]*10 - count_cnds_C[cdn_f])*(freq_codonC[cdn_nonf]/(1-freq_codonC[cdn_f]))		
			chiCr = (codon_cnt_droppedC[cdn_f,cdn_nonf] - Ecodon_cnt_droppedC[cdn_f,cdn_nonf])*(codon_cnt_droppedC[cdn_f,cdn_nonf] - Ecodon_cnt_droppedC[cdn_f,cdn_nonf])/Ecodon_cnt_droppedC[cdn_f,cdn_nonf]
			#file_objectC.write(f"{cdn_f},{cdn_nonf},{codon_cnt_droppedC[cdn_f,cdn_nonf]},{Ecodon_cnt_droppedC[cdn_f,cdn_nonf]},{chi}\n")
			chi_valsC.append(chiCr)
		
	chi_sumC = sum(chi_valsC)
	
	file_object_rand.write(f"{i},{chi_sumC}\n")
file_object_rand.close()	
								