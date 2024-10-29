This archive performs a comparison between codon enrichment in 5' v core for 1355 bacteria species and, both codon enrichment in highly vs lowly expressed genes as reported by Cambray et al, and those by Goodman et al

Cambray_logodds_HvL5.csv reports log odds ratios for each codon to show enrichment in highly vs lowly expressed genes from the Cambray et al dataset.
See directory TableS1 for full information on the preprocessing. It contains the same information as file Cambray_LO_PNInorm.csv

Reference sequence data was retrieved from NCBI RefSeq on 13th September 2023 (https://www.ncbi.nlm.nih.gov/refseq/).
See bact_accessions.csv for the full list of species analysed and their relative accession number.
Each .fna file containing the reference genome of one bacteria species per genus is not uploaded due to assembled file size but all data was then be placed in directory all_genomes to follow on with the analysis.


getSeq.py
	Scans the CDS of all genes to keep clean entries only and limit to one entry per gene name.	
	For each good quality gene: retrieves gene name, protein ID, number of exons.
	Also extracts 11 codons-long sequences from 5' end (including ATG) and gene core (around the middle of CDS).
	5' and core sequences output to two separate csv files in the directory "all_seqs".
	Script also calculates GC and GC3 content of each genes CDS from codon 11 onwards and outputs to separate file.
	Done for all bacteria files in all_genomes

enrichment.py
	Counts enrichment for each codon across all genes at 5' CDS and CDS core.
	Calculates log odds for each focal codon, comparing it to its synonyms (those coding for the same amino acid)
	Finds log odds ratio to compare 5' and core enrichment (log odds 5' CDS / log odds CDS core).
	For each genus, output to a single csv file with codon, amino acid, nucleotide at third position, log odds ratio (5' vs core) and standard error (SEM).
	Outputs placed in the directory "enrichment".
	Done for all bacteria files in all_seqs

merge_csv.r
	Combines within one file the Cambray log odds and all the log odds ratios for each codon from each bacterial genus.
	Output: all_log_odds.csv.

merge_csv_G.r
	Combines within one file the Goodman log odds and all the log odds ratios for each codon from each bacterial genus.
	Output: all_log_odds_G.csv.

getGCmean.R
	Finds mean GC and GC3 content across all genes in each genus.
	Output to single file: allbacteria_GC.csv

plotFig_bactcors.r
	Finds the Pearson correlation between the Goodman/Cambray log odds, and the log odds of any bacterial species.
	Produces a csv and histogram of the Pearson correlation values and the respective significance (with and without Bonferroni correction)
	Also plots the Pearson correlations against the average GC3 count for each genera
	Outputs: Fig9.pdf, correlation_df_G.csv (Goodman) and correlation_df_C.csv (Cambray) merged to TableS3.csv


