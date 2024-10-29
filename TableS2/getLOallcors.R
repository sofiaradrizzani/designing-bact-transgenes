install.packages("tableHTML", repos = "http://cran.us.r-project.org")
library(tableHTML)
library(Hmisc)


data1 <- read.csv("LO_stab_prot.csv", header = TRUE)
attach(data1)
df1 <- data.frame(data1)

df_LOstab <- subset(df1, select = -c(codon, amino_acid,ends_with))



#create a cor matrix with stars for significance
corstars <-function(x, method=c("pearson", "spearman"), result=c("none", "html")){
  #Compute correlation matrix
  require(Hmisc)
  require(tableHTML)
  
  x <- as.matrix(x)
  correlation_matrix<-rcorr(x, type=method[1])
  R <- correlation_matrix$r # Matrix of correlation coefficients
  p <- correlation_matrix$P # Matrix of p-value
  
  ## define notions for significance levels; spacing is important.
  mystars <- ifelse(p < .0001, "****", ifelse(p < .001, "*** ", ifelse(p < .01, "**  ", ifelse(p < .05, "*   ", "    "))))
  
  ## truncate the correlation matrix to two decimal
  R <- format(round(cbind(rep(-1.11, ncol(x)), R), 2))[,-1]
  
  ## build a new matrix that includes the correlations with their appropriate stars
  Rnew <- matrix(paste(R, mystars, sep=""), ncol=ncol(x))
  diag(Rnew) <- paste(diag(R), " ", sep="")
  rownames(Rnew) <- colnames(x)
  colnames(Rnew) <- paste(colnames(x), "", sep="")
  
  ## return the correlation matrix
  if (result[1]=="none") return(Rnew)
  if (result[1]=="html") return(tableHTML(Rnew))
}


cors_table <- corstars(df_LOstab, method="spearman", result="html")
#write_tableHTML(cors_table, file = "LOstabprot_cormatrix.html")
write_tableHTML(cors_table, file = "TableS2.html")

