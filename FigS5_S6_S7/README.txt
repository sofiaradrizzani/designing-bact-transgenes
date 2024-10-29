

Input data:
Cambray_DataS15.csv is the supplementary Data 15 directly from Cambray et al. Due to the excessive file size, this couldn't be uploaded. Please download (from https://osf.io/pm2qa?view_only=0d5b05fb08d84b76b21f399e832808b6) and rename file "Cambray_DataS15.csv" to progress with the analysis.
1241934tables1.csv is the supplementary directly from Goodman et al


split5seq.py
	Takes Cambray_DataS15.csv and 1241934tables1.csv as input to split the 5’ sequence into individual codons. 
	It then creates Cambray_split.csv and Goodman_split.csv containing all relevant information from the original files and, in addition, 11 new columns for each of the 11 codons in the 5’ sequence.

CalcLOHvL5_C.py and CalcLOHvL5_G.py
	Calculates log odds for each codon position for enrichment in PNI/RNAss (for Cambray) and Trans (for Goodman).
	Run for both Cambray_split.csv (_C) and Goodman_split.csv (_G)
	A new file is created for each codon position in LO_split directory (i.e. logodds_HvL_cod2.csv, …, logodds_HvL_cod11.csv)

mergeLObypos.R
	For both data sources, takes log odds ratios by position (i.e. logodds_HvL_cod2.csv, …, logodds_HvL_cod11.csv) and merges into a single file.
	Outputs: Goodman_LObypos.csv and Cambray_LObypos.csv and All_LObypos.csv (for both datasets)

plotFigS5.R
	Uses Goodman_LObypos.csv and Cambray_LObypos.csv created above to make 10 plots comparing Cambray v Goodman codon enrichment
	Output: FigS5.pdf

plotFigS6_heatmap.R 
	Uses Goodman_LObypos.csv and Cambray_LObypos.csv to generate a heat map showing where the log odds preferences match between the two data sources
	Output: FigS6.pdf

plotFigS7_bypos_tables.R
	Uses Goodman_LObypos.csv and Cambray_LObypos.csv to generate a pdf with two tables (for Goodman and Cambray) with the codon with highest log odds by 5' position
	Output: FigS7.pdf