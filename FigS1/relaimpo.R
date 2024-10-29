library(relaimpo)
install.packages("svglite", repos = "http://cran.us.r-project.org")
library(svglite)
library(gridExtra)

infile_C <- "GC_by_codon_Cambray.csv"
infile_G <- "GC_by_codon_Goodman.csv"

#Cambray	
df <- read.csv(infile_C)
	
colnames(df)[3:12] <- paste(c('p'), c(2:11), sep='_')
df <- df[,1:12]
	
#print(paste(fname, ":", dim(df)))

formula <- protein_level ~ p_2 + p_3 + p_4 +
                 p_5 + p_6 + p_7 + p_8 +
                 p_9 + p_10 + p_11

lm_model <- lm(formula, data = df)
	
relaimpo_results_C <- boot.relimp(lm_model, b = 1000, type = "lmg", rela=TRUE, rank = FALSE)
	
#print(booteval.relimp(relaimpo_results, sort = FALSE))
#pdf(outname, width = 4, height = 4)
	
	
#Goodman
df <- read.csv(infile_G)
	
colnames(df)[3:12] <- paste(c('p'), c(2:11), sep='_')
df <- df[,1:12]
	
#print(paste(fname, ":", dim(df)))

formula <- protein_level ~ p_2 + p_3 + p_4 +
                 p_5 + p_6 + p_7 + p_8 +
                 p_9 + p_10 + p_11

lm_model <- lm(formula, data = df)
	
relaimpo_results_G <- boot.relimp(lm_model, b = 1000, type = "lmg", rela=TRUE, rank = FALSE)


pdf("FigS1.pdf", width=12, height=10)

#pdf("FigS1.pdf", width=12, height=15)
par(mfrow = c(2,1))

plot(booteval.relimp(relaimpo_results_C, sort = FALSE, level=0.9), main ="A.")
plot(booteval.relimp(relaimpo_results_G, sort = FALSE, level=0.9), main ="B.")

dev.off()

