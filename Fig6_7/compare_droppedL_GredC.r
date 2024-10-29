library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
library(ggrepel)
library(ggpubr)

#read in VedIO

data1 = read.csv("Cambray_Goodman_codondrop_LOR.csv", h = TRUE)
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
	
fig6a = ggplot(df_dropG, aes(x = name, y = value)) + ylim(-0.1, 0.1) +
	geom_col(aes(fill = factor(ends_with, level = nts))) +  ggtitle("A.") +
	#scale_fill_manual(values = c("#EFE350FF", "#593D9CFF", "#B8627DFF", "#F68F46FF")) + 
	scale_fill_manual(values = c("#F68224", "#EFDD08","#53318E", "#4DA9AB")) +   
	geom_errorbar(aes(ymin = value - sem, ymax = value + sem, width = 0.4), colour = "black") + 
	facet_grid(~amino_acid, scales = "free_x", space = "free_x") + 
	scale_x_discrete(guide = guide_axis(angle = -90)) + 
	labs(x = "codon", y = "mean(log odds ratio without codon)", fill = "3rd Position \nNucleotide") + 
	theme(strip.text.x = element_text(face = "bold"), strip.background = element_rect(color = "black"),
	axis.title.x=element_text(size=11), axis.title.y=element_text(size=11),
axis.text.x=element_text(size=8), axis.text.y=element_text(size=8),
plot.title=element_text(size=11))
	
	
	
data1 = read.csv("Cambray_codondrop_LOR.csv", h = TRUE)
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
	
fig6b = ggplot(df_dropC, aes(x = name, y = value)) + ylim(-0.1, 0.1) +
	geom_col(aes(fill = factor(ends_with, level = nts))) +  ggtitle("B.") +
	#scale_fill_manual(values = c("#EFE350FF", "#593D9CFF", "#B8627DFF", "#F68F46FF")) + 
	scale_fill_manual(values = c("#F68224", "#EFDD08","#53318E", "#4DA9AB")) +   
	geom_errorbar(aes(ymin = value - sem, ymax = value + sem, width = 0.4), colour = "black") + 
	facet_grid(~amino_acid, scales = "free_x", space = "free_x") + 
	scale_x_discrete(guide = guide_axis(angle = -90)) + 
	labs(x = "codon", y = "mean(log odds ratio without codon)", fill = "3rd Position \nNucleotide") + 
	theme(strip.text.x = element_text(face = "bold"), strip.background = element_rect(color = "black"),
	axis.title.x=element_text(size=11), axis.title.y=element_text(size=11),
axis.text.x=element_text(size=8), axis.text.y=element_text(size=8),
plot.title=element_text(size=11))	
	



data2 <- read.csv("Goodman_logodds_high_low.csv", header = T)


attach(data2)

G_lor <- data2$log_odds
codonsG <- data2$codon

data3 <- read.csv("Cambray_logodds_high_low.csv", header = T)


attach(data3)

C_lor <- data3$log_odds
codonsC <- data3$codon

reg <- lm(G_lor~C_lor)
res <- reg$residuals


xl <- expression("Goodman predicted by Cambray, residuals")
yl <- expression("mean(log odds ratio without codon):Goodman")
df.xyres <- data.frame(y=res, x=df_dropG$value, codon = codonsC)
p.xyres <- ggplot(df.xyres, aes(x=x,y=y, label = codon)) + ggtitle("C.") +
  geom_point(size=0.8, color = "#53318E")+ theme(aspect.ratio=1) + 
  geom_smooth(method='lm', formula= y~x, colour = "#4DA9AB") +
  stat_cor(label.y = 1.75, label.x=-0.07, method = "spearman") +
  theme_bw()  + geom_text_repel(size=2) +labs(x = yl, y=xl)

df.xy <- data.frame(y=res, x=df_dropC$value, codon = codonsC)
xl <- expression("Goodman predicted by Cambray, residuals")
yl <- expression("mean(log odds ratio without codon):Cambray")

p.xy <- ggplot(df.xy, aes(x=x,y=y, label = codon)) + ggtitle("D.") +
  geom_point(size=0.8, color = "#53318E")+ theme(aspect.ratio=1) + 
  geom_smooth(method='lm', formula= y~x, colour = "#4DA9AB") +
  stat_cor(label.y = 1.75, label.x=0, method = "spearman") +theme_bw()  + 
  geom_text_repel(size=2) +labs(x = yl, y=xl)


pdf("Fig7.pdf", width =10, height = 15)
grid.arrange(
	grobs = list(fig6a, fig6b,  p.xyres, p.xy),
	widths = c(1,1),
	heights = c(1,1,1),
	layout_matrix = rbind(c(1,1), c(2,2), c(3,4)),
	respect=TRUE
	)
dev.off()