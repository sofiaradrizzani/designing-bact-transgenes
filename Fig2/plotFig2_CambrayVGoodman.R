
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
library(ggrepel)
library(ggpubr)
library(glue)



#read in log odds ratio
data_h = read.csv("TableS1.csv", h = TRUE)
attach(data_h)


#plot log odds ratio by codon, coloured based on 3rd position nucleotide
#Cambray
LO_C <- data_h$LO_PNInorm
class_cdn <- data_h$codon
aa <- data_h$amino_acid
ends_with <- data_h$ends_with

nts = c("A", "T", "G", "C")


df_C <- data.frame(
	amino_acid = aa,
	codon = toupper(class_cdn),
	valueC = LO_C,
	ends_with = ends_with
	)
	
fig2a <- ggplot(df_C, aes(x = codon, y = valueC))  +
	geom_col(aes(fill = factor(ends_with, level = nts))) +  ggtitle("A.") +
	#scale_fill_manual(values = c("#EFE350FF", "#593D9CFF", "#B8627DFF", "#F68F46FF")) +  
	scale_fill_manual(values = c("#F68224", "#EFDD08","#53318E", "#4DA9AB")) +   
	facet_grid(~amino_acid, scales = "free_x", space = "free_x") + 
	scale_x_discrete(guide = guide_axis(angle = -90)) + 
	labs(x = "codon", y = "log odds ratio", fill = "3rd Position \nNucleotide") + 
	theme(strip.text.x = element_text(face = "bold"), strip.background = element_rect(color = "black"))
	
#Goodman
LO_G <- data_h$VedIO

df_G <- data.frame(
	amino_acid = aa,
	codon = toupper(class_cdn),
	valueG = LO_G,
	ends_with = ends_with
	)

fig2b <- ggplot(df_G, aes(x = codon, y = valueG))  +
	geom_col(aes(fill = factor(ends_with, level = nts))) +  ggtitle("B.") +
	#scale_fill_manual(values = c("#EFE350FF", "#593D9CFF", "#B8627DFF", "#F68F46FF")) +  
	scale_fill_manual(values = c("#F68224", "#EFDD08","#53318E", "#4DA9AB")) +   
	facet_grid(~amino_acid, scales = "free_x", space = "free_x") + 
	scale_x_discrete(guide = guide_axis(angle = -90)) + 
	labs(x = "codon", y = "log odds ratio", fill = "3rd Position \nNucleotide") + 
	theme(strip.text.x = element_text(face = "bold"), strip.background = element_rect(color = "black"))
	


#XY plots
aa <- data_h$amino_acid

codon_pairs <- c()
amino_acid_pair <- c()
d.prots <- c()
d.exps <- c()

N <- length(aa)-1
y <- data_h$LO_PNInorm
x <- data_h$VedIO

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
  ggtitle("D.") +
  geom_point(size=0.8, color = "#53318E")+ 
  #theme(aspect.ratio=1, 
        #panel.background = element_blank(),
        #panel.grid.major = element_blank(), 
        #panel.grid.minor = element_blank(),
        #panel.border = element_rect(colour = "black", fill=NA)) + 
  geom_abline(slope=slp, intercept =int, col ="#F68224") +
  stat_cor(label.y = -1, label.x=2, method = "spearman") +
  theme_bw()  + 
  geom_text_repel(max.overlaps=nrow(df.xy),size=1.8) +
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
  ggtitle("C.") +
  geom_point(size=0.8, color = "#53318E")+ 
  #theme(aspect.ratio=1, 
        #panel.background = element_blank(),
        #panel.grid.major = element_blank(), 
        #panel.grid.minor = element_blank(),
        #panel.border = element_rect(colour = "black", fill=NA)) + 
  geom_abline(slope=slp, intercept =int, col ="#F68224") +
  stat_cor(label.y = -0.8, label.x=1, method = "spearman") +
  theme_bw()  + 
  geom_text_repel(max.overlaps= nrow(df.xy),size=1.8) +
  labs(x = xl, y=yl)
  

pdf("Fig2.pdf", width =10, height = 15)
grid.arrange(
	grobs = list(fig2a, fig2b, p.xy, p.diff),
	widths = c(1,1),
	heights = c(1,1,1),
	layout_matrix = rbind(c(1,1), c(2,2), c(3,4)),
	respect=TRUE
	)
dev.off()




