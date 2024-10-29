
library(glue)


#read in log odds ratio for each codon position
data_name <- glue("LO_split/Cambray/logodds_HvL5_cod2.csv")
data_cod2_C = read.csv(data_name, h = TRUE)
attach(data_cod2_C)
data_name <- glue("LO_split/Goodman/logodds_HvL5_cod2.csv")
data_cod2_G = read.csv(data_name, h = TRUE)
attach(data_cod2_G)

data_name <- glue("LO_split/Cambray/logodds_HvL5_cod3.csv")
data_cod3_C = read.csv(data_name, h = TRUE)
attach(data_cod3_C)
data_name <- glue("LO_split/Goodman/logodds_HvL5_cod3.csv")
data_cod3_G = read.csv(data_name, h = TRUE)
attach(data_cod3_G)

data_name <- glue("LO_split/Cambray/logodds_HvL5_cod4.csv")
data_cod4_C = read.csv(data_name, h = TRUE)
attach(data_cod4_C)
data_name <- glue("LO_split/Goodman/logodds_HvL5_cod4.csv")
data_cod4_G = read.csv(data_name, h = TRUE)
attach(data_cod4_G)

data_name <- glue("LO_split/Cambray/logodds_HvL5_cod5.csv")
data_cod5_C = read.csv(data_name, h = TRUE)
attach(data_cod5_C)
data_name <- glue("LO_split/Goodman/logodds_HvL5_cod5.csv")
data_cod5_G = read.csv(data_name, h = TRUE)
attach(data_cod5_G)

data_name <- glue("LO_split/Cambray/logodds_HvL5_cod5.csv")
data_cod5_C = read.csv(data_name, h = TRUE)
attach(data_cod5_C)
data_name <- glue("LO_split/Goodman/logodds_HvL5_cod5.csv")
data_cod5_G = read.csv(data_name, h = TRUE)
attach(data_cod5_G)

data_name <- glue("LO_split/Cambray/logodds_HvL5_cod6.csv")
data_cod6_C = read.csv(data_name, h = TRUE)
attach(data_cod6_C)
data_name <- glue("LO_split/Goodman/logodds_HvL5_cod6.csv")
data_cod6_G = read.csv(data_name, h = TRUE)
attach(data_cod6_G)

data_name <- glue("LO_split/Cambray/logodds_HvL5_cod7.csv")
data_cod7_C = read.csv(data_name, h = TRUE)
attach(data_cod7_C)
data_name <- glue("LO_split/Goodman/logodds_HvL5_cod7.csv")
data_cod7_G = read.csv(data_name, h = TRUE)
attach(data_cod7_G)

data_name <- glue("LO_split/Cambray/logodds_HvL5_cod8.csv")
data_cod8_C = read.csv(data_name, h = TRUE)
attach(data_cod8_C)
data_name <- glue("LO_split/Goodman/logodds_HvL5_cod8.csv")
data_cod8_G = read.csv(data_name, h = TRUE)
attach(data_cod8_G)

data_name <- glue("LO_split/Cambray/logodds_HvL5_cod9.csv")
data_cod9_C = read.csv(data_name, h = TRUE)
attach(data_cod9_C)
data_name <- glue("LO_split/Goodman/logodds_HvL5_cod9.csv")
data_cod9_G = read.csv(data_name, h = TRUE)
attach(data_cod9_G)

data_name <- glue("LO_split/Cambray/logodds_HvL5_cod10.csv")
data_cod10_C = read.csv(data_name, h = TRUE)
attach(data_cod10_C)
data_name <- glue("LO_split/Goodman/logodds_HvL5_cod10.csv")
data_cod10_G = read.csv(data_name, h = TRUE)
attach(data_cod10_G)

data_name <- glue("LO_split/Cambray/logodds_HvL5_cod11.csv")
data_cod11_C = read.csv(data_name, h = TRUE)
attach(data_cod11_C)
data_name <- glue("LO_split/Goodman/logodds_HvL5_cod11.csv")
data_cod11_G = read.csv(data_name, h = TRUE)
attach(data_cod11_G)




#merge into single file and export csv for each data source
df_LObypos_C <- data.frame(data_cod2_C["codon"],data_cod2_C["amino_acid"],
	data_cod2_C["log_odds"], data_cod3_C["log_odds"],data_cod4_C["log_odds"],data_cod5_C["log_odds"],
	data_cod6_C["log_odds"],data_cod7_C["log_odds"],data_cod8_C["log_odds"],data_cod9_C["log_odds"],
	data_cod10_C["log_odds"],data_cod11_C["log_odds"])
	
colnames(df_LObypos_C) <- c("codon","amino_acid","2","3","4","5","6","7","8","9","10","11")

outfile_name <- "Cambray_LObypos.csv"
write.csv(df_LObypos_C, file=outfile_name, row.names=FALSE)


df_LObypos_G <- data.frame(data_cod2_G["codon"],data_cod2_G["amino_acid"],
	data_cod2_G["log_odds"], data_cod3_G["log_odds"],data_cod4_G["log_odds"],data_cod5_G["log_odds"],
	data_cod6_G["log_odds"],data_cod7_G["log_odds"],data_cod8_G["log_odds"],data_cod9_G["log_odds"],
	data_cod10_G["log_odds"],data_cod11_G["log_odds"])
	
colnames(df_LObypos_G) <- c("codon", "amino_acid","2","3","4","5","6","7","8","9","10","11")

outfile_name <- "Goodman_LObypos.csv"
write.csv(df_LObypos_G, file=outfile_name, row.names=FALSE)


df_both <- data.frame(data_cod2_C["codon"],data_cod2_C["amino_acid"],
						data_cod2_G["log_odds"],data_cod2_C["log_odds"],
						data_cod3_G["log_odds"],data_cod3_C["log_odds"],
						data_cod4_G["log_odds"],data_cod4_C["log_odds"],
						data_cod5_G["log_odds"],data_cod5_C["log_odds"],
						data_cod6_G["log_odds"],data_cod6_C["log_odds"],
						data_cod7_G["log_odds"],data_cod7_C["log_odds"],
						data_cod8_G["log_odds"],data_cod8_C["log_odds"],
						data_cod9_G["log_odds"],data_cod9_C["log_odds"],
						data_cod10_G["log_odds"],data_cod10_C["log_odds"],
						data_cod11_G["log_odds"],data_cod11_C["log_odds"])

colnames(df_both) <- c("codon", "amino_acid","G2", "C2","G3", "C3","G4", "C4","G5", "C5",
						"G6", "C6","G7", "C7","G8", "C8","G9", "C9","G10", "C10","G11", "C11")

outfile_name <- "All_LObypos.csv"
write.csv(df_both, file=outfile_name, row.names=FALSE)
