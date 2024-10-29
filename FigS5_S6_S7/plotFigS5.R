
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
library(ggrepel)
library(ggpubr)
library(glue)


#read in dataframes (already aligned)
cdf <- read.csv('Cambray_LObypos.csv', header=TRUE)
colnames(cdf)[3:12] <- c(2:11)

gdf <- read.csv('Goodman_LObypos.csv', header=TRUE)
colnames(gdf)[3:12] <- c(2:11)



#initialise empty dataframe
df_all <- data.frame(matrix(ncol=4,nrow=0, dimnames=list(NULL, c("codon", "cam_logodds", "good_logodds", "position"))))

#create dataframe in long form
for(n in 3:12){
  cam_logodds <- cdf[,n]
  good_logodds <- gdf[,n]
  codon <- cdf[,1]
  ccol_pos <- colnames(cdf[n])
  ccol_name <- paste(c('position'), ccol_pos, sep='')
  temp_df <- data.frame(cbind(codon, cam_logodds, good_logodds))
  temp_df$position <- ccol_name
  df_all <- data.frame(rbind(df_all, temp_df))
}

df_all <- na.omit(df_all)
df_all$cam_logodds <- as.numeric(df_all$cam_logodds)
df_all$good_logodds <- as.numeric(df_all$good_logodds)
df_all$codon <- toupper(df_all$codon)



#define function that will generate plots for each given position
#uses method 2 pca for line, otherwise strange
plotter <- function(n){
  pos_name <- paste(c('position'), n, sep='')
  pos_title <- paste(c('Position'), n, sep=' ')
  df_sec <- df_all[df_all$position == pos_name,]
  df_sec$status <- c('normal')
  df_sec[df_sec$codon == 'AGG',]$status <- 'funky'
  cols <- c('normal' = '#53318E', 'funky'='#4DA9AB')
  
  pca_res <- prcomp(df_sec[, c('cam_logodds', 'good_logodds')], scale.=TRUE)
  slope <- pca_res$rotation[2, 1] / pca_res$rotation[1, 1]
  mean_x <- mean(df_sec[['cam_logodds']], na.rm = TRUE)
  mean_y <- mean(df_sec[['good_logodds']], na.rm = TRUE)
  intercept <- mean_y - slope * mean_x
  res <- ggplot(df_sec, aes(x=cam_logodds, y=good_logodds, label=codon)) + 
    #geom_point(aes(color=status), pch=21) + 
    geom_point(aes(color=status)) + 
    theme_bw() + 
    labs(title=paste(pos_title), x=c('Cambray codon enrichment'), y=c('\nGoodman codon enrichment')) + 
    geom_abline(slope = slope, intercept = intercept, color = "#F68224", linewidth=0.7) +  
    stat_cor(method="spearman", label.x.npc = "left") + geom_text_repel(size=3) + 
    scale_color_manual(values = cols) + 
    theme(legend.position='none', 
          axis.title.x=element_text(size=18), axis.title.y=element_text(size=18),
          axis.text.x=element_text(size=13), axis.text.y=element_text(size=13),
          plot.title=element_text(size=18))
  
  return(res)
}


#empty list of plots
plotz <- list()
#fill empty list of plots
for(n in 2:11){
  plowt_name <- paste(c('Position'), n, sep='')
  assign(plowt_name, plotter(n))
  plotz[[plowt_name]] <- plotter(n)
}

#put out plot list to pdf
pdf("FigS5.pdf", width=15, height=35)
grid.arrange(grobs = plotz, ncol=2, nrow=5)
dev.off()


