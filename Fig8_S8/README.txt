
This archive performed a dropped log odds analysis but based on RNAfold stability measures rather than protein/RNA

Input data:
Cambray_DataS15.csv is the supplementary Data 15 directly from Cambray et al. Due to the excessive file size, this couldn't be uploaded. Please download (from https://osf.io/pm2qa?view_only=0d5b05fb08d84b76b21f399e832808b6) and rename file "Cambray_DataS15.csv" to progress with the analysis.
1241934tables1.csv is the supplementary directly from Goodman et al


STABILITY_CALCULATOR.py
	Requires input files Cambray_DataS15.csv and 1241934tables1.csv. 
	Filters the Cambray data (and outputs to extracted_data.csv)
	Uses the package ViennaRNA (specifically the RNAfold function) to calculate the mRNA stability (MFE) of the 5' region (2-11) of the sequences from the Goodman and Cambray data.
	Outputs: GSTR.csv (Goodman) and CSTR.csv (Cambray).

Calc_LO_RNAfold_C.py
	Finds the log odds ratios in the high vs low stability RNAfold measures (from CSTR.csv)
	It also finds the average log odds of high vs low stability when each 5' codon is dropped (i.e. the average log odds of the other 9 positions in the 10-codon 5' sequence)
	Outputs: Cambray_LOstab.csv and Cambray_LOstab_drop.csv

Calc_LO_RNAfold_G.py
	Same as Calc_LO_RNAfold_C.py but for Goodman et al data
	Outputs: Goodman_LOstab.csv and Goodman_LOstab_drop.csv

G_C_drop_analysis.py
	Finds the log odds with dropped codons for the sequences in the Goodman data (1241934tables1.csv) but using the Cambray_LOstab.csv rather than the log odds calculated from the Goodman data.
	Outputs: Cambray_Goodman_LOstab_drop.csv

plotFigS8_LOdrop_GwithC.R
	Panel A: bar plot of the Cambray derived log odds ratios applied to the Goodman constructs (Cambray_Goodman_LOstab_drop.csv)
	Panel B: bar plot of the Cambray derived log odds ratios applied to the Cambray constructs (Cambray_LOstab_drop.csv)
	Panel C: XY plot of Cambray_LOstab_drop.csv v Goodman_LOstab_drop.csv
	Panel D: XY plot of the Cambray and Goodman log odds ratios, with regression to derive the residuals
	Panel E: plot of the calculated residuals by the mean of mean dropped log odds
	Output: FigS8.pdf

plotFig8_LOdrop_GwithG.R
	Panel A: bar plot of the Goodman derived log odds ratios applied to the Goodman constructs (Goodman_LOstab_drop.csv)
	Panel B: bar plot of the Cambray derived log odds ratios applied to the Cambray constructs (Cambray_LOstab_drop.csv)
	Panel C: XY plot of Cambray_LOstab_drop.csv v Goodman_LOstab_drop.csv
	Panel D: XY plot of the Cambray and Goodman log odds ratios, with regression to derive the residuals
	Panel E: plot of the calculated residuals by the mean of mean dropped log odds
	Output: Fig8.pdf
