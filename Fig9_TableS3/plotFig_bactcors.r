library(dplyr)
library(ggplot2)
library(broom)
library(glue)



###Goodman
#read in csv with VPNInorm + V5core for all 1355 bacteria
log_odds_data = read.csv("all_log_odds_G.csv", h = TRUE)

#extract column with Cambray log odds
exp.lo = log_odds_data$VedIO
ncols = dim(log_odds_data)[2]
r.vals = c()
p.vals = c()

#calculate r and p values for each genome's V5-core correlated against VedIO
for (i in c(5:ncols)) {
	lo = unlist(log_odds_data[i])
	ct = cor.test(exp.lo, lo, method = "pearson")
	
	ct.r = ct$estimate
	ct.r = round(ct.r, digits=4)
	
	ct.p = ct$p.value
	ct.p = round(ct.p, digits=4)

	r.vals = c(r.vals, ct.r)
	p.vals = c(p.vals, ct.p)
}

#find the points that demarcate statistical significance 
accession = colnames(log_odds_data)[5:ncols]
df.cor = data.frame(acc = accession, rvals = r.vals, pvals = p.vals)
df.cor_G = data.frame(species = accession, Goodman_r = r.vals, Goodman_pvalue = p.vals) #just for csv output purposes
non_sig = which(df.cor$pvals > 0.05)
df.non_sig = df.cor$rvals[non_sig]
max.r = max(df.non_sig)
min.r = min(df.non_sig)

non_sig_bonf = which(df.cor$pvals > (0.05/ncols)) #includes bonferroni correction for multiple testing
df.non_sig_bonf = df.cor$rvals[non_sig_bonf]
max.r_bonf = max(df.non_sig_bonf)
min.r_bonf = min(df.non_sig_bonf)



#make output csv
write.csv(x <- df.cor, file = "correlation_df_G.csv")

#merge GC data and Pearson correlation data
corr_data <- read.csv("correlation_df_G.csv", h = TRUE)
gc_data <- read.csv("allbacteria_GC.csv", h = TRUE)

merged_data <- merge(corr_data, gc_data, by = "acc")
merged_data2 <- transform(merged_data, GC3_squared=GC3_mean^2)

quadratic_model <- lm(merged_data2$rvals ~ merged_data2$GC3_mean + merged_data2$GC3_squared)
sum_model <- glance(quadratic_model)
adj_r <- sum_model$adj.r.squared
adj_r <- round(adj_r,digits=3)
pval <- sum_model$p.value
pval <- signif(pval, digits=3)


#make figure

#plot histogram of R values 
ghist <- ggplot(corr_data, aes(x=rvals)) + 
geom_histogram(binwidth=0.02, fill = "#53318E") +
xlab("Pearson's r") +
ylab("Frequency") +
ggtitle("A.")

phist_G <- ghist + 
xlim(-0.50,0.85) +
ylim(0,150) +
theme(aspect.ratio=1, 
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA)) +
  geom_vline(xintercept = max.r, col = "#F68224") +
  geom_vline(xintercept = min.r, col = "#F68224") +
  geom_vline(xintercept = max.r_bonf, col = "#4DA9AB") +
  geom_vline(xintercept = min.r_bonf, col = "#4DA9AB") +
  geom_vline(xintercept = 0.7553, col = "black", lty=2) #E.coli r=0.7553, p=0


#plot GC3 v Pearson correlation
pGC3_G <- ggplot(merged_data, aes(x = GC3_mean, y = rvals)) +
	geom_point(colour = "#53318E", pch = 21, size = 0.5) + 
	theme(aspect.ratio=1, 
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA)) +
	#geom_point(data = ecoli, colour = "#F2824B", pch = 18, cex = 4) +
	#geom_text(data = ecoli, label = "E. coli", nudge_x = 5, nudge_y = 0.05, colour = "#F2824B") +
	#theme_bw() + 
	geom_smooth(method = "lm", formula = y ~ x + I(x^2), colour = "#F68224", se = FALSE, linewidth = 0.7) +
	annotate("text", x=78, y=-0.1, label= glue("Adjusted R-squared: {adj_r},\np-value: {pval}"), size = unit(3, "pt")) +
	labs(x = "GC3 content (%)", y = "Pearson's r") +
	ggtitle("B.")
	
	
	
	

###Cambray
#read in csv with VPNInorm + V5core for all 1355 bacteria
log_odds_data = read.csv("all_log_odds.csv", h = TRUE)

#extract column with Cambray log odds
exp.lo = log_odds_data$VPNInorm
ncols = dim(log_odds_data)[2]
r.vals = c()
p.vals = c()

#calculate r and p values for each genome's V5-core correlated against VPNInorm
for (i in c(5:ncols)) {
	lo = unlist(log_odds_data[i])
	ct = cor.test(exp.lo, lo, method = "pearson")
	
	ct.r = ct$estimate
	ct.r = round(ct.r, digits=4)
	
	ct.p = ct$p.value
	ct.p = round(ct.p, digits=4)

	r.vals = c(r.vals, ct.r)
	p.vals = c(p.vals, ct.p)
}

#find the points that demarcate statistical significance 
accession = colnames(log_odds_data)[5:ncols]
df.cor = data.frame(acc = accession, rvals = r.vals, pvals = p.vals)
df.cor_C = data.frame(species = accession, Cambray_r = r.vals, Cambray_pvalue = p.vals) #just for csv output purposes
non_sig = which(df.cor$pvals > 0.05)
df.non_sig = df.cor$rvals[non_sig]
max.r = max(df.non_sig)
min.r = min(df.non_sig)

non_sig_bonf = which(df.cor$pvals > (0.05/ncols)) #includes bonferroni correction for multiple testing
df.non_sig_bonf = df.cor$rvals[non_sig_bonf]
max.r_bonf = max(df.non_sig_bonf)
min.r_bonf = min(df.non_sig_bonf)



#make output csv
write.csv(x <- df.cor, file = "correlation_df_C.csv")

#merge GC data and Pearson correlation data
corr_data <- read.csv("correlation_df_C.csv", h = TRUE)
gc_data <- read.csv("allbacteria_GC.csv", h = TRUE)

merged_data <- merge(corr_data, gc_data, by = "acc")
merged_data2 <- transform(merged_data, GC3_squared=GC3_mean^2)

quadratic_model <- lm(merged_data2$rvals ~ merged_data2$GC3_mean + merged_data2$GC3_squared)
sum_model <- glance(quadratic_model)
adj_r <- sum_model$adj.r.squared
adj_r <- round(adj_r,digits=3)
pval <- sum_model$p.value
pval <- signif(pval, digits=3)


#make figure

#plot histogram of R values 
ghist <- ggplot(corr_data, aes(x=rvals)) + 
geom_histogram(binwidth=0.02, fill = "#53318E") +
xlab("Pearson's r") +
ylab("Frequency") +
ggtitle("C.")

phist_C <- ghist + 
xlim(-0.50,0.85) +
ylim(0,150) +
theme(aspect.ratio=1, 
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA)) +
  geom_vline(xintercept = max.r, col = "#F68224") +
  geom_vline(xintercept = min.r, col = "#F68224") +
  geom_vline(xintercept = max.r_bonf, col = "#4DA9AB") +
  geom_vline(xintercept = min.r_bonf, col = "#4DA9AB") +
  geom_vline(xintercept = 0.504, col = "black", lty=2) #E.coli r=0.504, p=0


#plot GC3 v Pearson correlation
pGC3_C <- ggplot(merged_data, aes(x = GC3_mean, y = rvals)) +
	geom_point(colour = "#53318E", pch = 21, size = 0.5) + 
	theme(aspect.ratio=1, 
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA)) +
	#geom_point(data = ecoli, colour = "#F2824B", pch = 18, cex = 4) +
	#geom_text(data = ecoli, label = "E. coli", nudge_x = 5, nudge_y = 0.05, colour = "#F2824B") +
	#theme_bw() + 
	geom_smooth(method = "lm", formula = y ~ x + I(x^2), colour = "#F68224", se = FALSE, linewidth = 0.7) +
	annotate("text", x=78, y=-0.1, label= glue("Adjusted R-squared: {adj_r},\np-value: {pval}"), size = unit(3, "pt")) +
	labs(x = "GC3 content (%)", y = "Pearson's r") +
	ggtitle("D.")




pdf("Fig9.pdf", height = 10, width = 10)

gridExtra::grid.arrange(phist_G, pGC3_G, phist_C, pGC3_C, 
                        ncol = 2, nrow = 2)
                        
dev.off()



#add accession column then output merged correlation data to csv
merged_corr_data <- merge(df.cor_G, df.cor_C, by = "species")

bact_acc <- read.csv("bact_accessions.csv", h = TRUE)
corr_acc_data <- merge(bact_acc, merged_corr_data, by = "species")

write.csv(x <- corr_acc_data, file = "TableS3.csv")


#what proportion of genomes show significant correlation?
#tot_non_sig = length(non_sig)
#tot = dim(df.cor)[1]
#prop_sig = (tot - tot_non_sig) / tot

#which genomes show significant + correlation and which significant - correlation?
#sig_pos = df.cor$acc[df.cor$rvals > 0 & df.cor$pvals < (0.05/ncols)]
#sig_neg = df.cor$acc[df.cor$rvals < 0 & df.cor$pvals < (0.05/ncols)]



