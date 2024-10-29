
This archive performed a dropped log odds analysis based on protein/RNA measures

Input data:
Cambray_DataS15.csv is the supplementary Data 15 directly from Cambray et al. Due to the excessive file size, this couldn't be uploaded. Please download (from https://osf.io/pm2qa?view_only=0d5b05fb08d84b76b21f399e832808b6) and rename file "Cambray_DataS15.csv" to progress with the analysis.
1241934tables1.csv is the supplementary directly from Goodman et al

Cambray_analysis.py and Goodman_analysis.py
	Finds the log odds ratios in the high vs low PNI/RNAss (for Cambray) and ProtFCC (for Goodman)
	It also finds the average log odds of high vs low measures when each 5' codon is dropped (i.e. the average log odds of the other positions in the 10-codon 5' sequence)
	Outputs: Cambray_logodds_high_low.csv, Cambray_codondrop_LOR.csv, Cambray_codondrop_LOR.csv and Goodman_logodds_high_low.csv, Goodman_codondrop_LOR.csv, Goodman_codondrop_LOR.csv
	
compare_droppedL.r
	Takes Goodman_codondrop_LOR.csv and Cambray_codondrop_LOR.csv to plot the log odds ratios derived from each respective dataset
	Takes Goodman_logodds_high_low.csv and Cambray_logodds_high_low.csv to consider the regression between the two. It then finds the residuals from their regression line and sees whether these are predicted by the mean log odds of the flanking codons for both the Goodman data and the Cambray data.
	Output: Fig6.pdf
	
Cambray_Good_drop_analysis.py
	Same as above but uses 1241934tables1.csv and Cambray_logodds_high_low.csv so that for Goodman constructs the log odds ratios derived from the Cambray dataset are applied rather than the ones from Goodman
	Output: Cambray_Goodman_codondrop_LOR.csv

compare_droppedL_GredC.r
	Same as compare_droppedL.r but uses Cambray_Goodman_codondrop_LOR.csv rather than Goodman_codondrop_LOR.csv
	Output: Fig7.pdf
