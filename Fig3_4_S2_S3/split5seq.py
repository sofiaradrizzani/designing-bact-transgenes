import pandas as pd
import re
import os

cwd = os.getcwd()

#create output files
outfile_path_C = os.path.join(cwd, "Cambray_split.csv")
outfile_C = open(outfile_path_C, "w")
outfile_C.write("id,PNI,PNInorm,seq,cod2,cod3,cod4,cod5,cod6,cod7,cod8,cod9,cod10,cod11\n")

outfile_path_G = os.path.join(cwd, "Goodman_split.csv")
outfile_G = open(outfile_path_G, "w")
#outfile_G.write("id,protFCC,trans,seq,cod2,cod3,cod4,cod5,cod6,cod7,cod8,cod9,cod10,cod11\n")
outfile_G.write("id,protFCC,seq,cod2,cod3,cod4,cod5,cod6,cod7,cod8,cod9,cod10,cod11\n")


#create dataframes from csv, and then to list
#Cambray
columns = ["id","gs.sequence", "clean.lin.prot.mean", "ss.rna.dna.mean"]
dfC = pd.read_csv("Cambray_DataS15.csv", usecols = columns)
dfC = dfC.rename(columns={"id":"id" ,"clean.lin.prot.mean": "clean_lin_prot_mean","gs.sequence": "gs_sequence", "ss.rna.dna.mean": "ss_rna_dna_mean"})
dfC['clean_lin_prot_mean_norm'] = dfC['clean_lin_prot_mean']/dfC['ss_rna_dna_mean'] #find protein by rna
dfC = dfC.dropna()
dfC_as_list = dfC.values.tolist()

#Goodman
columns = ["Name","CDS.seq", "Prot.FCC"]
#columns = ["Name","CDS.seq", "Prot.FCC", "Trans"]
dfG = pd.read_csv("1241934tables1.csv", usecols = columns)
dfG = dfG.rename(columns={"Name":"Name","Prot.FCC": "ProtFCC","CDS.seq": "CDSseq", "Trans": "Trans"})
dfG = dfG.rename(columns={"Name":"Name","Prot.FCC": "ProtFCC","CDS.seq": "CDSseq"})
dfG = dfG.dropna()
dfG_as_list = dfG.values.tolist()


#Cambray: iterate through each row and get required info (including split 5' seq) -> write out to new csv
for row in dfC_as_list:
    id = row[0]
    five_seq = row[1]
    PNI = row[2]
    PNInorm = row[4]

    five_seq = five_seq.upper()

    #print(five_seq)
    n = 3
    if re.search('[A-Za-z]{11}', str(five_seq)):
        split_seq = [five_seq[i:i + n] for i in range(0, 30, n)]


    #cod1 = split_seq[0]
    cod2 = split_seq[0]
    cod3 = split_seq[1]
    cod4 = split_seq[2]
    cod5 = split_seq[3]
    cod6 = split_seq[4]
    cod7 = split_seq[5]
    cod8 = split_seq[6]
    cod9 = split_seq[7]
    cod10 = split_seq[8]
    cod11 = split_seq[9]

    #print(cod1)
    #print(cod2)

    outfile_C.write(f"{id},{PNI},{PNInorm},{five_seq},{cod2},{cod3},{cod4},{cod5},{cod6},{cod7},{cod8},{cod9},{cod10},{cod11}\n")

outfile_C.close()

#Goodman: iterate through each row and get required info (including split 5' seq) -> write out to new csv
for row in dfG_as_list:
    id = row[0]
    five_seq = row[1]
    protFCC = row[2]
    #prot = row[3]
    #trans = row[2]

    # print(five_seq)
    n = 3
    if re.search('[A-Za-z]{11}', str(five_seq)):
        split_seq = [five_seq[i:i + n] for i in range(0, 33, n)]

    cod1 = split_seq[0]
    cod2 = split_seq[1]
    cod3 = split_seq[2]
    cod4 = split_seq[3]
    cod5 = split_seq[4]
    cod6 = split_seq[5]
    cod7 = split_seq[6]
    cod8 = split_seq[7]
    cod9 = split_seq[8]
    cod10 = split_seq[9]
    cod11 = split_seq[10]

    # print(cod1)
    # print(cod2)

    #outfile_G.write(f"{id},{prot},{trans},{five_seq},{cod2},{cod3},{cod4},{cod5},{cod6},{cod7},{cod8},{cod9},{cod10},{cod11}\n")
    outfile_G.write(f"{id},{protFCC},{five_seq},{cod2},{cod3},{cod4},{cod5},{cod6},{cod7},{cod8},{cod9},{cod10},{cod11}\n")

outfile_G.close()