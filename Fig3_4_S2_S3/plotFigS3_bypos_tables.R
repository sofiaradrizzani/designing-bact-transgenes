library(wrapr)
library(tidyverse)
library("gridExtra")
library(ggplot2)
library(ggpubr)
library(hash)
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


aminoacids <- unique(dfc$Amino_acid)
alphabet <- toupper(letters)



cmat <- matrix(, nrow = length(aminoacids), ncol = 10)

for(n in 2:11){
position_name <- paste(c('position'), n, sep='')
df_sec <- dfc[dfc$position == position_name,]
xindex <- n-1
for(aa in aminoacids){
df_secc <- df_sec[df_sec$Amino_acid == aa,]
cmox <- df_secc[df_secc$cam_logodds == max(df_secc$cam_logodds), 1]
yindex <- match(c(aa), aminoacids)
cmat[yindex, xindex] <- c(cmox)
}
}


camdf <- as.data.frame(as.table(cmat))
colnames(camdf) <- c("Amino_acid", "position", "top_codon")
camdf$Amino_acid <- match(camdf$Amino_acid, alphabet)
camdf$Amino_acid <- aminoacids[camdf$Amino_acid]
camdf$position <- match(camdf$position, alphabet)
camdf$position <- as.numeric(camdf$position)
camdf$position <- camdf$position + 1

cova1 <- ggplot(camdf, aes(x = factor(position), y = factor(Amino_acid))) +
  geom_tile(color='black', fill='white') +
  geom_text(aes(label = top_codon), size = 4) +
  scale_fill_manual(values = c("Y" = "white", "N" = "black")) +
  labs(title = 'Cambray', x = "Position", y = "Amino Acid") +
  theme_minimal() +
  theme(legend.position = "none")




gmat <- matrix(, nrow = length(aminoacids), ncol = 10)

for(n in 2:11){
position_name <- paste(c('position'), n, sep='')
df_sec <- dfc[dfc$position == position_name,]
xindex <- n-1
for(aa in aminoacids){
df_secc <- df_sec[df_sec$Amino_acid == aa,]
gmox <- df_secc[df_secc$good_logodds == max(df_secc$good_logodds), 1]
yindex <- match(c(aa), aminoacids)
if(length(gmox) == 1){
gmat[yindex, xindex] <- c(gmox)
} else {
gmat[yindex, xindex] <- c('NA')
}
}
}



gooddf <- as.data.frame(as.table(gmat))
colnames(gooddf) <- c("Amino_acid", "position", "top_codon")
gooddf$Amino_acid <- match(gooddf$Amino_acid, alphabet)
gooddf$Amino_acid <- aminoacids[gooddf$Amino_acid]
gooddf$position <- match(gooddf$position, alphabet)
gooddf$position <- as.numeric(gooddf$position)
gooddf$position <- gooddf$position + 1

gova1 <- ggplot(gooddf, aes(x = factor(position), y = factor(Amino_acid))) +
  geom_tile(color='black', fill='white') +
  geom_text(aes(label = top_codon), size = 4) +
  scale_fill_manual(values = c("Y" = "white", "N" = "black")) +
  labs(title = 'Goodman', x = "Position", y = "Amino Acid") +
  theme_minimal() +
  theme(legend.position = "none")



#pdf('FIGURE_4_SUPPL_PNI_RNAss.pdf', width=10, height=20)
pdf('FigS3.pdf', width=10, height=20)
grid.arrange(gova1, cova1, ncol=1, nrow=2)
dev.off()

