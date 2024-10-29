library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
library(ggrepel)
library(ggpubr)




data1 = read.csv("Cambray_native_LO_drop.csv", h = TRUE)
attach(data1)

lor <- data1$mean_LOR_not
class_cdn <- data1$cdn
sem <- data1$sem
sd <- data1$sd
aa <- data1$aa

nts = c("A", "T", "G", "C")

df_dropG <- data.frame(
	amino_acid = aa,
	name = class_cdn,
	value = lor,
	sem = sem
	)
	
figB = ggplot(df_dropG, aes(x = name, y = value)) + ylim(-0.1, 0.1) +
	geom_col(aes(fill = factor(ends_with, level = nts))) +  
	ggtitle("B.") +
	#scale_fill_manual(values = c("#EFE350FF", "#593D9CFF", "#B8627DFF", "#F68F46FF")) + 
	scale_fill_manual(values = c("#F68224", "#EFDD08","#53318E", "#4DA9AB")) +   
	geom_errorbar(aes(ymin = value - sem, ymax = value + sem, width = 0.4), colour = "black") + 
	facet_grid(~amino_acid, scales = "free_x", space = "free_x") + 
	scale_x_discrete(guide = guide_axis(angle = -90)) + 
	labs(x = "codon", y = "mean(log odds ratio without codon)", 
	fill = "3rd Position \nNucleotide") + 
	theme(strip.text.x = element_text(face = "bold"), strip.background = element_rect(color = "black"),
	axis.title.x=element_text(size=12), axis.title.y=element_text(size=12),
	plot.title=element_text(size=12))
	
	
	
data1 = read.csv("Cambray_LO_drop.csv", h = TRUE)
attach(data1)

lor <- data1$mean_LOR_not
class_cdn <- data1$cdn
sem <- data1$sem
sd <- data1$sd
aa <- data1$aa


df_dropC <- data.frame(
	amino_acid = aa,
	name = class_cdn,
	value = lor,
	sem = sem
	)

figA = ggplot(df_dropC, aes(x = name, y = value)) + ylim(-0.1, 0.1) +
	geom_col(aes(fill = factor(ends_with, level = nts))) +  
	ggtitle("A.") +
	#scale_fill_manual(values = c("#EFE350FF", "#593D9CFF", "#B8627DFF", "#F68F46FF")) + 
	scale_fill_manual(values = c("#F68224", "#EFDD08","#53318E", "#4DA9AB")) +   
	geom_errorbar(aes(ymin = value - sem, ymax = value + sem, width = 0.4), colour = "black") + 
	facet_grid(~amino_acid, scales = "free_x", space = "free_x") + 
	scale_x_discrete(guide = guide_axis(angle = -90)) + 
	labs(x = "codon", y = "mean(log odds ratio without codon)", 
	fill = "3rd Position \nNucleotide") + 
	theme(strip.text.x = element_text(face = "bold"), strip.background = element_rect(color = "black"),
	axis.title.x=element_text(size=12), axis.title.y=element_text(size=12),
	plot.title=element_text(size=12))	
	



data2 <- read.csv("native_LO_drop.csv", h = TRUE)
attach(data2)


data3 <- read.csv("Cambray_LO_drop.csv", h = TRUE)
attach(data3)

codonsC <- data3$cdn


df.xy <- data.frame(x=data2$mean_LOR_not, y=data3$mean_LOR_not, codon = codonsC)
xl <- expression("native log odds ratios")
yl <- expression("Cambray log odds ratios")

p.xy <- ggplot(df.xy, aes(x=x,y=y, label = codon)) + 
ggtitle("C.") +
  geom_point(size=0.8, color = "#53318E")+ theme(aspect.ratio=1) + 
  theme(aspect.ratio=1,
  axis.title.x=element_text(size=12), axis.title.y=element_text(size=12),
  plot.title=element_text(size=12))+ 
  geom_smooth(method='lm', formula= y~x, colour= "#4DA9AB") +
  stat_cor(label.y = 0.06, label.x=-0.015, method = "spearman") +
  theme_bw()  + 
  geom_text_repel(max.overlaps= nrow(df.xy), size=2) +
  labs(x = yl, y=xl)


pdf("FigS9.pdf", width =10, height = 15)
grid.arrange(
	grobs = list(figA, figB, p.xy),
	widths = c(1,1),
	heights = c(1,1,1),
	layout_matrix = rbind(c(1,1), c(2,2), c(3,NA)),
	respect=TRUE
	)
dev.off()





