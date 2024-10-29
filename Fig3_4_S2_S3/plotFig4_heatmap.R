library(wrapr)
library(tidyverse)
library("gridExtra")
library(ggplot2)
library(ggpubr)
install.packages("hash", repos = "http://cran.us.r-project.org")
library(hash)
install.packages("naniar", repos = "http://cran.us.r-project.org")
library(naniar)


#Read in csvs with positional logodds vectors for goodman and cambray data
cdf <- read.csv('Cambray_LObypos.csv', header=TRUE)
colnames(cdf)[3:12] <- c(2:11)

gdf <- read.csv('Goodman_LObypos.csv', header=TRUE)
colnames(gdf)[3:12] <- c(2:11)


#initialise empty dataframe
dfc <- data.frame(matrix(ncol=5,nrow=0, dimnames=list(NULL, c("codon", "amino_acid", "cam_logodds", "good_logodds", "position"))))

#get dataframe in long format
for(n in 2:11){
cam_logodds <- cdf[,n+1]
good_logodds <- gdf[,n+1]
codon <- gdf[,1]
Amino_acid <- gdf[,2]
ccol_pos <- colnames(cdf[n+1])
ccol_name <- paste(c('position'), ccol_pos, sep='')
temp_df <- data.frame(cbind(codon, Amino_acid, cam_logodds, good_logodds))
temp_df$position <- ccol_name
dfc <- data.frame(rbind(dfc, temp_df))
}


dfc$cam_logodds <- as.numeric(dfc$cam_logodds)
dfc$good_logodds <- as.numeric(dfc$good_logodds)
dfc$codon <- toupper(dfc$codon)

#create matrix of Y and N (agree and disagree on codon with highest log-odds ratio in synonymous block) where x axis is position 2-11 and y axis is a list of amino acids.

aminoacids <- unique(dfc$Amino_acid)

mat <- matrix(, nrow = length(aminoacids), ncol = 10)


for(n in 2:11){
position_name <- paste(c('position'), n, sep='')
df_sec <- dfc[dfc$position == position_name,]
xindex <- n-1
for(aa in aminoacids){
df_secc <- df_sec[df_sec$Amino_acid == aa,]
cmox <- df_secc[df_secc$cam_logodds == max(df_secc$cam_logodds), 1]
gmox <- df_secc[df_secc$good_logodds == max(df_secc$good_logodds), 1]
yindex <- match(c(aa), aminoacids)
if(length(gmox) == 1){
if(cmox == gmox){
mat[yindex, xindex] <- c('Y')
}
if(cmox != gmox){
mat[yindex, xindex] <- c('N')
}
} else {
mat[yindex, xindex] <- NA
}
}
}

#create long df from matrix to plot out
alphabet <- toupper(letters)

df <- as.data.frame(as.table(mat))
colnames(df) <- c("Amino_acid", "position", "value")
df$Amino_acid <- match(df$Amino_acid, alphabet)
df$Amino_acid <- aminoacids[df$Amino_acid]
df$position <- match(df$position, alphabet)
df$position <- as.numeric(df$position)
df$position <- df$position + 1




ova <- ggplot(df, aes(x = factor(position), y = factor(Amino_acid), fill=value)) +
  geom_tile(color = "black") +
  #scale_fill_manual(values = c("Y" = "#F68224", "N" = "#53318E"), na.value="gray") +
  #scale_fill_manual(values = c("Y" = "white", "N" = "#53318E"), na.value="gray") +
  scale_fill_manual(values = c("Y" = "white", "N" = "#212121"), na.value="gray") +
  labs(x = "Position", y = "Amino Acid") +
  #labs(title='Agreement', x = "Position", y = "Amino Acid") +
  geom_text(aes(label = ifelse(is.na(value), "NA", '')),
            color = "black", size = 4) +
  theme_minimal() +
  theme(
    legend.position = "none")



#output to pdf
pdf('Fig4.pdf', width=10, height=10)
grid.arrange(ova, ncol=1, nrow=1)
dev.off()

#############################################################################
#Get data for binomial tests for each amino acid row


syn_codons <- list(
Phe = c('ttt', 'ttc'),
Ser = c('tct', 'tcc', 'tca', 'tcg', 'agt', 'agc'),
Tyr = c('tat', 'tac'),
Cys = c('tgt', 'tgc'),
Leu = c('ctt', 'ctc', 'cta', 'ctg', 'tta', 'ttg'),
Pro = c('cct', 'ccc', 'cca', 'ccg'),
His = c('cat', 'cac'),
Gln = c('caa', 'cag'),
Arg = c('cgt', 'cgc', 'cga', 'cgg', 'aga', 'agg'),
Ile = c('att', 'atc', 'ata'),
Thr = c('act', 'acc', 'aca', 'acg'),
Asn = c('aat', 'aac'),
Lys = c('aaa', 'aag'),
Val = c('gtt', 'gtc', 'gta', 'gtg'),
Ala = c('gct', 'gcc', 'gca', 'gcg'),
Asp = c('gat', 'gac'),
Glu = c('gaa', 'gag'),
Gly = c('ggt', 'ggc', 'gga', 'ggg')
)


three_to_one_letter <- c(
  Ala = "A", Arg = "R", Asn = "N", Asp = "D",
  Cys = "C", Gln = "Q", Glu = "E", Gly = "G",
  His = "H", Ile = "I", Leu = "L", Lys = "K",
  Met = "M", Phe = "F", Pro = "P", Ser = "S",
  Thr = "T", Trp = "W", Tyr = "Y", Val = "V"
)


amino_acids <- sapply(names(syn_codons), function(x) three_to_one_letter[[x]])

statdf <- data.frame(Amino_acid = amino_acids, degen = sapply(syn_codons, length))

statdf$agrees <- 'h'


for(amino in amino_acids){
df_am <- df[df$Amino_acid == amino,]
df_am <- na.omit(df_am)
agreez <-  sum(df_am$value == 'Y')
statdf[statdf$Amino_acid == amino,]$agrees <- agreez
}

statdf$null_conc <- 1/statdf$degen
statdf$agrees <- as.numeric(statdf$agrees)


statdf$nas <- 'n'

for(amino in amino_acids){
df_am <- df[df$Amino_acid == amino,]
#nas <-  sum(df_am$value == NA)
nas <- sum(is.na(df_am$value))
statdf[statdf$Amino_acid == amino,]$nas <- nas
}


pwals <- c()

for(n in 1:nrow(statdf)){
row <- statdf[n,]
aa <- df[df$Amino_acid == row$Amino_acid,]
aam <- na.omit(aa)
naam <- nrow(aam)
X <- row$agrees
N <- naam
p <- row$null_conc
pwal <- binom.test(X, N, p)$p.value
pwals <- c(pwals, pwal)
}


pwows <- p.adjust(pwals, method = "bonferroni")

statdf$pvalsBon <- pwows
statdf$pvals <- pwals

statdf$nas <- as.numeric(statdf$nas)
statdf$weight <- ((10 - statdf$nas) / 10)

statdf$inds <- (statdf$null_conc * statdf$weight) / sum(statdf$weight)

weight_mean_null_conc <- sum(statdf$inds)

#mean_null_conc <- mean(statdf$null_conc)
overall_t <- binom.test(sum(statdf$agrees), length(na.omit(df$value)), weight_mean_null_conc)


#for X^2
#statdf$exp <- statdf$null_conc * (nrow(statdf)*10)
#chi.test()
