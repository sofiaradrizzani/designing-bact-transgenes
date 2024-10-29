library(wrapr)
library(tidyverse)
library("gridExtra")
library(ggplot2)
library(ggpubr)
library(ggrepel)


#read in dataframe
df_all <- read.csv('All_LObypos.csv', header=TRUE)
attach(df_all)


plotz <- list()

for(n in 2:11){
plot_name <- paste(c('plot'), n, sep='')
plot_title <- paste(c('Position'), n, sep=' ')

#find LO of the position for each dataset -> make new dataframe
df_sel <- data.frame(df_all$codon, df_all$amino_acid)
colnames(df_sel) <- c("codon", "amino_acid")
LO_G <- df_all[,n+(n-1)] 
LO_C <- df_all[,n+n]
df_sel$good_LO <- LO_G
df_sel$camb_LO <- LO_C

#df_sel <- na.omit(df_sel)

#find codon pairs and calculate log odds differences
aa <- df_sel$amino_acid

codon_pairs <- c()
amino_acid_pair <- c()
d.prots <- c()
d.exps <- c()

N <- length(aa)-1
y <- df_sel$good_LO
x <- df_sel$camb_LO

df.xy <- data.frame(x=x, y=y)

for (i in c(1:N)) {
		st <- i + 1
		nd <- length(aa)
	
	for (j in c(st:nd)) {
		prot1 <- y[i]
		exp1 <- x[i]
		aa1 <- aa[i]
		cdn1 <- codon[i]
		
		prot2 <- y[j]
		exp2 <- x[j]
		aa2 <- aa[j]
		cdn2 <- codon[j]
		
		#print(c(i, j, aa1, aa2, cdn1, cdn2))
			
		if (identical(aa1,aa2)==TRUE) {
		
			if (exp1 > exp2) {
			d.exp <- exp1 - exp2
			cdn_pair <- paste(cdn1,cdn2, sep =":")
			d.prot <- prot1 -prot2
			} else {
			d.exp <- exp2 - exp1
			cdn_pair <- paste(cdn2,cdn1, sep =":")
			d.prot <- prot2 -prot1
			}
			
			codon_pairs <- c(codon_pairs, cdn_pair)
			amino_acid_pair <- c(amino_acid_pair, aa1)
			d.prots <- c(d.prots, d.prot)
			d.exps <- c(d.exps, d.exp)
			}

}
}

df.diff <- data.frame(amino_acid_p = amino_acid_pair, codon_pair = codon_pairs, deltaLO_C = d.prots, deltaLO_G = d.exps)

attach(df.diff)


pca <- prcomp(~deltaLO_C+deltaLO_G, df.diff)
slp <- with(pca, rotation[2,1] / rotation[1,1])
int <- with(pca, center[2] - slp*center[1])
yl <- expression(Delta~Cambray~codon~enrichment)
xl <- expression(Delta~Goodman~codon~enrichment)
options(ggrepel.max.overlaps = Inf)

p.diff <- ggplot(df.diff, aes(x=deltaLO_C,y=deltaLO_G, label = codon_pair)) +   
  ggtitle(plot_title) +
  geom_point(size=0.8, color = "#53318E")+ 
  #theme(aspect.ratio=1, 
        #panel.background = element_blank(),
        #panel.grid.major = element_blank(), 
        #panel.grid.minor = element_blank(),
        #panel.border = element_rect(colour = "black", fill=NA)) + 
  geom_abline(slope=slp, intercept =int, col ="#F68224") +
  stat_cor(method = "spearman") +
  #stat_cor(label.y = -1.25, label.x=1, method = "spearman") +
  theme_bw()  + 
  geom_text_repel(max.overlaps= nrow(df.xy),size=1.8) +
  labs(x = xl, y=yl)

#plot_name <- p.diff
#plotz <- c(plotz,plot_name)

assign(plot_name, p.diff)
plotz[[plot_name]] <- p.diff

detach(df.diff)

}

#export to pdf
pdf("FigS2.pdf", width=15, height=35)
grid.arrange(grobs = plotz, ncol=2, nrow=5)
dev.off()







