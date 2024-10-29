This archive contains scripts and files necessary to analyse passenger effects (non-random distributions of codons) amongst the Cambray and Goodman data


Directory getNulls:
Input data: Cambray_logodds_high_low.csv and Goodman_logodds_high_low.csv are the log odds ratios of enrichment in highly vs lowly expressed transgenes.
getNull.py takes the input data and generates a null of random association for focal codon and any given non focal codon. Outputs: Cambray_null.csv and Goodman_null.csv
The sd of the null Cambray vector is manually transferred to make_hist_2plots.R


Directory sd_test:
Input data:
	Cambray_logodds_high_low.csv (as above)
	Cambray_DataS15.csv is the supplementary Data 15 directly from Cambray et al. Due to the excessive file size, this couldn't be uploaded. Please download (from https://osf.io/pm2qa?view_only=0d5b05fb08d84b76b21f399e832808b6) and rename file "Cambray_DataS15.csv" to progress with the analysis.
	1241934tables1.csv is the supplementary Table 1 directly from Goodman et al
	Cambray_Goodman_codondrop_LOR.csv contains dropped log odds scores (see directory Fig6_7 for further details)
Cambray_Mantel_analysis.py uses this input data to find the standard deviations between the mean log odds of the non-focal codons per construct. Output: Cambray_codondrop_LOR_Mantel_sds.csv


Directory nonrandom_analysis:
Nonrandom_analysis2.py takes 1241934tables1.csv and Cambray_DataS15.csv. Outputs: Goodman_chi.csv and Rand_chi_Efromfreq.csv


Directory make_fig:
make_hist_2plots.R
	Takes Cambray_codondrop_LOR_Mantel_sds.csv and Rand_chi_Efromfreq.cav from above, as well as Cambray_codondrop_LOR.csv (see directory Fig6_7 for further details).
	Output:Fig5.pdf

