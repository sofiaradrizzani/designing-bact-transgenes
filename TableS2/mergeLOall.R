

#read in Cambray stability LO data files
data_C_PNInorm = read.csv("Cambray_LO_PNInorm.csv", h = TRUE)
attach(data_C_PNInorm)

data_C_utr = read.csv("Cambray_LOstab_utr.csv", h = TRUE)
attach(data_C_utr)

data_C_fivep = read.csv("Cambray_LOstab_fivep.csv", h = TRUE)
attach(data_C_fivep)

data_C_RNAfold = read.csv("Cambray_LOstab_RNAfold.csv", h = TRUE)
attach(data_C_RNAfold)


#read in Goodman stability LO data files
data_G_prot = read.csv("Goodman_LO_protFCC.csv", h = TRUE)
attach(data_G_prot)

data_G_trans = read.csv("Goodman_LO_trans.csv", h = TRUE)
attach(data_G_trans)

data_G_dG = read.csv("Goodman_LOstab_dG.csv", h = TRUE)
attach(data_G_dG)

data_G_dGnoutr = read.csv("Goodman_LOstab_dGnoutr.csv", h = TRUE)
attach(data_G_dGnoutr)

data_G_dGunif = read.csv("Goodman_LOstab_dGunif.csv", h = TRUE)
attach(data_G_dGunif)

data_G_RNAfold = read.csv("Goodman_LOstab_RNAfold.csv", h = TRUE)
attach(data_G_RNAfold)


#merge into single file and export to single csv
df_LOstab <- data.frame(data_C_PNInorm["codon"], data_C_PNInorm["amino_acid"], data_C_PNInorm["ends_with"],
	data_C_PNInorm["LO_PNInorm"], 
	data_C_utr["LO_utrCdsStructureMFE"], data_C_fivep["LO_fivepCdsStructureMFE"], data_C_RNAfold["LO_RNAfold_C"],
	data_G_prot["LO_protFCC"], data_G_trans["LO_trans"], 
	data_G_dG["LO_dG"], data_G_dGnoutr["LO_dGnoutr"], data_G_dGunif["LO_dGunif"], data_G_RNAfold["LO_RNAfold_G"])
	
write.csv(df_LOstab, file= "LO_stab_prot.csv", row.names=FALSE)

