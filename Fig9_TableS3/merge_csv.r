#script purpose: merge the expression log odds data with the log odds
#data for each genome into one dataframe 

expression_data = read.csv("Cambray_logodds_HvL5.csv", h = TRUE)

log_odds_df = expression_data

log_odds_df <- with(log_odds_df, log_odds_df[order(amino_acid, codon) , ])


path1 = file.path(getwd(), "enrichment", "*.csv")

files_list = Sys.glob(path1)

for (fl in files_list) {
	genome_data = read.csv(fl, h = TRUE)
	genome_data <- with(genome_data, genome_data[order(amino_acid, codon) , ])

	log_odds = genome_data[3]
	log_odds_df = cbind(log_odds_df, log_odds)
	#genome_data = data.frame(genome_data[1],genome_data[3])
	#merged_data = merge(expression_data, genome_data, by = "codon")

}


write.csv(x = log_odds_df, file = "all_log_odds.csv")

