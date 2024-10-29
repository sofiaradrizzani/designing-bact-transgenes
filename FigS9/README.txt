
This archive plots the equivalent of Fig6 but comparing log odds between Cambray and those based on RefSeq sequences and PaxDB protein measures (instead of Cambray v Goodman et al data)

Input data:
Cambray_DataS15.csv is the supplementary Data 15 directly from Cambray et al. Due to the excessive file size, this couldn't be uploaded. Please download (from https://osf.io/pm2qa?view_only=0d5b05fb08d84b76b21f399e832808b6) and rename file "Cambray_DataS15.csv" to progress with the analysis.
Ecoli.fna are sequences from the reference E coli genome obtained from NCBI RefSeq (accession: GCF_000008865.2 downloaded 13th September 2023)

getSeq.py
	Scans the CDS of all genes to keep clean entries only and limit to one entry per gene name.	
	For each good quality gene: retrieves gene name, protein ID, number of exons.
	Also extracts 11 codons-long sequences from 5' end (including ATG)
	Output: Ecoli_five.csv

Calc_LO_C.py
	Finds the log odds ratios in the high vs low expression taking protein/RNA measures (from Cambray_DataS15.csv)
	It also finds the average log odds of high vs low stability when each 5' codon is dropped (i.e. the average log odds of the other 9 positions in the 10-codon 5' sequence)
	Outputs: Cambray_LO.csv and Cambray_LO_drop.csv
	
Calc_LO_native.py
	For native sequences: finds the average log odds of high vs low stability when each 5' codon is dropped (i.e. the average log odds of the other 9 positions in the 10-codon 5' sequence)
	Outputs: native_LO_drop.csv

CambrayVnative_drop_analysis.py
	Finds the log odds with dropped codons for the sequences in the native sequences (NCBI RefSeq .csv) using Cambray log odds scores (Cambray_LOstab.csv).
	Outputs: Cambray_native_LO_drop.csv

plotFig_LOdrop_CambrayVnative.R
	Panel A: bar plot of the Cambray derived log odds ratios applied to the Cambray constructs (Cambray_LO_drop.csv)
	Panel B: bar plot of the Cambray derived log odds ratios applied to the native constructs (Cambray_native_LO_drop.csv)
	Panel C: plot of Cambray_LO_drop.csv v native_LO_drop.csv codon enrichment vectors
	Output: FigS9.pdf
	
	



