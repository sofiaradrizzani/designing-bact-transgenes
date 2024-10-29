library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
library(ggrepel)
library(ggpubr)


data2 <- read.csv("Cambray_codondrop_LOR.csv",header = T)
attach(data2)
vals <- unlist(data2$mean_LOR_not)
C_sd_fullset = sd(vals)


data1 <- read.csv("Cambray_codondrop_LOR_Mantel_sds.csv", header = T)

attach(data1)

sds <- unlist(data1$sd)

sds.obs <- sds[1]


#from python script getNull.py
null <- 0.006555275572412394

#now get all rands but observed
sds.rand <- sds[-1]

xmax <- max(sds.rand)
    

z = round((sds.obs - mean(sds.rand))/sd(sds.rand), 3)
textforplot = paste('Z = ', z, sep="")

z_null_C = round((null- mean(sds.rand))/sd(sds.rand), 3)
textforplot_null = paste('Z = ', z_null_C, sep="")


df <- data.frame(x = sds.rand)

p8b <- ggplot(df, aes(x=x)) + 
#geom_histogram(colour=4, fill= "#53318E", bins=50) + 
#geom_histogram(colour="#9060B0", fill= "#53318E", bins=50) + 
geom_histogram(colour="#53318E", fill= "#53318E", bins=50) + 
geom_vline(xintercept = sds.obs, color="#4DA9AB", linewidth=0.7) + 
ggtitle("B.")+ 
geom_vline(xintercept = C_sd_fullset, color="#EFDD08", linewidth=0.7) + 
geom_vline(xintercept = null, color="#F68224", linewidth=0.7) + 
theme(aspect.ratio=1) +theme_classic() + 
xlab("SD of codon dropout vectors") +ylab("Count")+ 
annotate("text", x = sds.obs-0.002, y=200, label = textforplot)+ 
annotate("text", x = null+0.002, y=200, label = textforplot_null)
  


data1a <- read.csv("Rand_chi_Efromfreq.csv", header = T)

attach(data1a)

chi <- unlist(data1a$chi)

chi.obs <- chi[1]

#  from theory E(chi) = df(chi), with 3660 
nullchi <- 3660 -1
#from chi squared table get 0.05 sign with 

one_percent = qchisq(p = .01, df = nullchi, lower.tail = FALSE)


#now get all rands but observed
chi.rand <- chi[-1]

xmax <- max(chi.rand)
    

z = round((chi.obs - mean(chi.rand))/sd(chi.rand), 3)
textforplotchi = paste('Z = ', z, sep="")

z_null_C = round((nullchi- mean(chi.rand))/sd(chi.rand), 3)
textforplot_nullchi = paste('Z = ', z_null_C, sep="")


df <- data.frame(x = chi.rand)

p8a <- ggplot(df, aes(x=x)) + 
#geom_histogram(colour=4, fill= "#53318E", bins=100) + 
#geom_histogram(colour="#9060B0", fill= "#53318E", bins=50) + 
geom_histogram(colour="#53318E", fill= "#53318E", bins=50) + 
geom_vline(xintercept = chi.obs, color="#4DA9AB", linewidth=0.7) + 
#geom_vline(xintercept = nullchi, color="#EFDD08")+ 
ggtitle("A.")+ 
geom_vline(xintercept = one_percent, color="#F68224", linewidth=0.7) + 
theme(aspect.ratio=1) +
theme_classic() + 
xlab("chi squared value") +ylab("Count")+ 
annotate("text", x = chi.obs-20000, y=700, label = textforplotchi)+ 
annotate("text", x = 22000, y=700, label = textforplot_nullchi)
  


pdf("Fig5.pdf", height = 5, width = 10)

grid.arrange(
	grobs = list(p8a, p8b),
	widths = c(1,1),
	heights = c(1),
	layout_matrix = rbind(c(1,2)),
	respect=TRUE
	)

dev.off()


