library(wrapr)
library(tidyverse)
library("gridExtra")
library(ggplot2)
library(ggpubr)
library(ggrepel)



#Read in input file
data_h <- read.csv('TableS1.csv', header=TRUE)

pdf("FigS4.pdf", width = 12, height = 6)

attach(data_h)

aa <- data_h$amino_acid

codon_pairs <- c()
amino_acid_pair <- c()
d.prots <- c()
d.exps <- c()

N <- length(aa)-1
y <- data_h$LO_PNInorm
x <- data_h$LO_trans

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

df.diff <- data.frame(amino_acid_p = amino_acid_pair, codon_pair = codon_pairs, delta_Vprot = d.prots, delta_VedIO = d.exps)

attach(df.diff)



pca <- prcomp(~delta_VedIO+delta_Vprot, df.diff)
slp <- with(pca, rotation[2,1] / rotation[1,1])
int <- with(pca, center[2] - slp*center[1])
#xl <- expression(Delta~V[Five_prime])
#yl <- expression(Delta~V[Core])
xl <- expression(Delta~Goodman~codon~enrichment)
yl <- expression(Delta~Cambray~codon~enrichment)
options(ggrepel.max.overlaps = Inf)

p.diff <- ggplot(df.diff, aes(x=delta_VedIO,y=delta_Vprot, label = codon_pair)) +   
  ggtitle("B.") +
  geom_point(size=0.8, color = "#53318E")+ 
  theme(aspect.ratio=1) +
  #theme(aspect.ratio=1,        
  		#panel.background = element_blank(),
        #panel.grid.major = element_blank(), 
        #panel.grid.minor = element_blank(),
        #panel.border = element_rect(colour = "black", fill=NA)) + 
  geom_abline(slope=slp, intercept =int, col ="#F68224") +
  #stat_cor(label.y = -1.25, label.x=1, method = "spearman") +
  stat_cor(method = "spearman") +
  theme_bw()  + 
  geom_text_repel(max.overlaps= nrow(df.xy),size=1.8) +
  labs(x = xl, y=yl)
  
pca <- prcomp(~x+y, df.xy)
slp <- with(pca, rotation[2,1] / rotation[1,1])
int <- with(pca, center[2] - slp*center[1])
#xl <- expression(V[Five_prime])
#yl <- expression(V[core])
xl <- "Goodman codon enrichment"
yl <- "Cambray codon enrichment"
options(ggrepel.max.overlaps = Inf)

p.xy <- ggplot(df.xy, aes(x=x,y=y, label = codon)) + 
  ggtitle("A.") +
  geom_point(size=0.8, color = "#53318E")+ 
  theme(aspect.ratio=1) +
  #theme(aspect.ratio=1, 
        #panel.background = element_blank(),
        #panel.grid.major = element_blank(), 
        #panel.grid.minor = element_blank(),
        #panel.border = element_rect(colour = "black", fill=NA)) + 
  geom_abline(slope=slp, intercept =int, col ="#F68224") +
  #stat_cor(label.y = -0.75, label.x=-1, method = "spearman") +
  stat_cor(method = "spearman") +
  theme_bw()  + 
  geom_text_repel(max.overlaps= nrow(df.xy),size=1.8) +
  labs(x = xl, y=yl)
  
	
grid.arrange(
	grobs = list(p.xy, p.diff),
	layout_matrix = rbind(c(1,2))
	)
dev.off()

