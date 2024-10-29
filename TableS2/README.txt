
This archive inspects stability measures and the compares the log odds derived by each of them

Input data:
Cambray_DataS15.csv is the supplementary Data 15 directly from Cambray et al. Due to the excessive file size, this couldn't be uploaded. Please download (from https://osf.io/pm2qa?view_only=0d5b05fb08d84b76b21f399e832808b6) and rename file "Cambray_DataS15.csv" to progress with the analysis.
1241934tables1.csv is the supplementary directly from Goodman et al


STABILITY_CALCULATOR.py
	Requires input files Cambray_DataS15.csv and 1241934tables1.csv. 
	Filters the Cambray data (and outputs to extracted_data.csv)
	Uses the package ViennaRNA (specifically the RNAfold function) to calculate the mRNA stability (MFE) of the 5' region (2-11) of the sequences from the Goodman and Cambray data.
	Outputs: GSTR.csv (Goodman) and CSTR.csv (Cambray).

Calc_allLO_C.py
	Calculates the log odds for the Cambray data for the RNAfold measures, and their columns gs.utrCdsStructureMFE, gs.fivepCdsStructureMFE and the division of columns clean.lin.prot.mean by ss.rna.dna.mean
	Output: Cambray_LOstab_UTRto30.csv, Cambray_LOstab_0to60.csv, Cambray_LOstab_RNAfold.csv, and Cambray_LO_PNInorm.csv

Calc_allLO_G.py
	Calculates the log odds for the Goodman data for the RNAfold measures, and their columns dG, dG.noutr, dG.unif and Prot.FCC, Trans
	Outputs: Goodman_LOstab_dG.csv, Goodman_LOstab_dGnoutr.csv, Goodman_LOstab_dGunif.csv, Goodman_LOstab_RNAfold.csv, and Goodman_LO_protFCC.csv, Goodman_LO_trans.csv

mergeLOall.R
	Merges all stability log odds from both data sources into a single file, along with the log odds calculated for protein or protein/RNA measures
	Outputs: LO_stab_prot.csv

getLOallcors.R
	Generates a correlation matrix comparing the log odds calculated from all stability measures between Cambray and Goodman data
	Significance is displayed with a ranging number of stars so that: p < .0001 = "****"", p < .001 = "*** ", p < .01 = "**  ", p < .05 = "*   ", p > .05 = no stars
	Output: TableS2.html



