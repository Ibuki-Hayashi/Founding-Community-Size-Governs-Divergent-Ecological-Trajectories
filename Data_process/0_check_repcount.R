#-- information of directory
fd <- "Output"

#-- read data from previous dir
source("Functions/functions.R")
seqtab <- readRDS("Data/merge_ASVtab.rds")
OTUtab <- readRDS("Data/seqOTUtab_99.rds")
taxtable <- readRDS("Data/merge_ASVtaxonomylist.rds")
sq99 <- read.table("Data/ASV_OTU_corestab_99.txt",header=TRUE)

exp_r <- read.csv("Data/expdata_rep.csv")
exp_i <- read.csv("Data/expdata_ino.csv")

#-- ASV 2 OTU by correspondence table
A2O <- function(ASVtable, cortab, key){
  cor <- cortab[cortab[,2] %in% key, 4:ncol(cortab)]
  rownames(cor) <- key
  u_cor <- cor[(rownames(cor) %in% colnames(ASVtable)),]
  u_cor <- u_cor[,(colSums(u_cor) != 0)]
  
  ans <- matrix(NA,ncol=ncol(u_cor),nrow=nrow(ASVtable))
  colnames(ans) <- colnames(u_cor); rownames(ans) <- rownames(ASVtable)
  for(i in 1:ncol(ans)){
    key <- rownames(u_cor)[(u_cor[,i] == 1)]
    if(length(key) == 1){
      ans[,i] <- ASVtable[,(colnames(ASVtable) %in% key)]
    }else{
      ans[,i] <- rowSums(ASVtable[,(colnames(ASVtable) %in% key)])
    }
  }
  return(ans)
}
OTUtab <- A2O(ASVtable=seqtab, cortab=sq99, key=rownames(taxtable)[taxtable[,3] != "c_"])
OTU_t <- taxtable[(1:ncol(OTUtab)),]

###### seqtab ###########
#- without stdDNA
stdos <- rownames(taxtable)[taxtable[,3] != "c_"]
rep_t <- taxtable[(rownames(taxtable) %in% stdos),]
rep_s <- cbind(Sample_name=rownames(seqtab[(str_sub(rownames(seqtab),1,5) == "HS3_S"), (colnames(seqtab) %in% stdos)]), as.data.frame(seqtab[(str_sub(rownames(seqtab),1,5) == "HS3_S"), (colnames(seqtab) %in% stdos)]))
rep_w <- cbind(Sample_name=rownames(seqtab[(str_sub(rownames(seqtab),1,5) == "HS3_W"), (colnames(seqtab) %in% stdos)]), as.data.frame(seqtab[(str_sub(rownames(seqtab),1,5) == "HS3_W"), (colnames(seqtab) %in% stdos)]))

rep_s_OTU <- cbind(Sample_name=rownames(OTUtab[(str_sub(rownames(OTUtab),1,5) == "HS3_S"),]), as.data.frame(OTUtab[(str_sub(rownames(OTUtab),1,5) == "HS3_S"),]))
rep_w_OTU <- cbind(Sample_name=rownames(OTUtab[(str_sub(rownames(OTUtab),1,5) == "HS3_W"),]), as.data.frame(OTUtab[(str_sub(rownames(OTUtab),1,5) == "HS3_W"),]))

D_soil <- merge(exp_r, rep_s, by="Sample_name"); D_water <- merge(exp_r, rep_w, by="Sample_name")
D_soil_OTU <- merge(exp_r, rep_s_OTU, by="Sample_name"); D_water_OTU <- merge(exp_r, rep_w_OTU, by="Sample_name")

#-- rarefaction curve
rr_func <- function(data, taxt, day, sna, sn, seed=1112){
  data <- data[(data$Day == day),]
  set.seed(seed); data <- data[sample(c(1:nrow(data)), size=sn),]
  du <- as.matrix(data[,(7:ncol(data))])

  pdf(sna, width=4, height=4)
  rarecurve(du, step=50, label=F)
  abline(v = 3000, col = "red", lty = 2, lwd = 2)
  dev.off()
}

if(length(list.dirs(sprintf("%s/FigS1", fd), recursive=F))==0){
  dir.create(sprintf("%s/FigS1", fd))
}

rr_func(data=D_soil, taxt=, day=2,
        sna=sprintf("%s/FigS1/Rarecurve_soil_ASV_day2.pdf", fd), sn=50)
rr_func(data=D_soil, taxt=, day=8,
        sna=sprintf("%s/FigS1/Rarecurve_soil_ASV_day8.pdf", fd), sn=50)

rr_func(data=D_water, taxt=, day=2,
        sna=sprintf("%s/FigS1/Rarecurve_water_ASV_day2.pdf", fd), sn=50)
rr_func(data=D_water, taxt=, day=8,
        sna=sprintf("%s/FigS1/Rarecurve_water_ASV_day8.pdf", fd), sn=50)

