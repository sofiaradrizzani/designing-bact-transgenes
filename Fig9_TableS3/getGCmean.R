library(glue)

home_dir <- getwd()

setwd("./all_seqs")
curr_dir <- getwd()

allfiles <- list.files(curr_dir, pattern = "_GC.csv", full.names = TRUE)


df_out <- data.frame(acc = character(), "GC_mean" = numeric(), "GC3_mean" = numeric())
                     
  for (f in allfiles) {
    genus <- sub(".*seqs/", "", f)
    genus <- sub("_GC.*", "", genus)
  
    data1 <- read.csv(f, header = TRUE)
    attach(data1)
    df1 <- data.frame(data1)
    
    
    mean_GC <- colMeans(data1["GC"])
    mean_GC3 <- colMeans(data1["GC3"])
    
    output <- c(genus,mean_GC,mean_GC3)
    
    df_out[nrow(df_out) + 1, ] <- output
      
    detach(data1)
    
    }
    

outfile_name <- glue("{home_dir}/allbacteria_GC.csv")

write.csv(df_out, outfile_name, row.names=FALSE)









