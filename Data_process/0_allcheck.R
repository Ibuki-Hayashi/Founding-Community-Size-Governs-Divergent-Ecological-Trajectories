#-- information of directory
fd <- "Output"

#-- read data from previous dir
source("Functions/functions.R")
seqtab <- readRDS("Data/merge_ASVtab.rds")
OTUtab <- readRDS("Data/seqOTUtab_99.rds")
taxtable <- readRDS("Data/merge_ASVtaxonomylist.rds")
sq99 <- read.table("Data/ASV_OTU_corestab_99.txt",header=TRUE)
sq97 <- read.table("Data/ASV_OTU_corestab_97.txt",header=TRUE)

exp_r <- read.csv("Data/expdata_rep.csv")
exp_i <- read.csv("Data/expdata_ino.csv")

#-- ASV 2 OTU by correspondence table
A2O <- function(ASVtable, cortab){
  key <- colnames(cortab[,4:ncol(cortab)])
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

###### seqtab ###########
#- rarefuction by n reads(without stdDNA)
raref <- function(seqtab,taxtable,n,seed=1111){
  std_out <- seqtab[,(colnames(seqtab) %in% rownames(taxtable[(str_sub(taxtable[,2],1,3) != "STD"),]))]
  use <- std_out[(rowSums(std_out)>(n-1)),]
  set.seed(seed)
  ans1 <- rrarefy(use,n); rownames(ans1) <- rownames(use)
  ans2 <- sprintf("ASV_seqtab_%s.rds",as.character(n))
  ans3 <- taxtable[(rownames(taxtable) %in% colnames(ans1)),]
  ans4 <- sprintf("taxonomylist_%s.rds",as.character(n))
  return(list(ans1,ans2,ans3,ans4))
}

number <- 3000 ## the threshold of rarefaction

seq_ans <- raref(seqtab=seqtab,taxtable=taxtable,n=number)

O97seq <- A2O(ASVtable=seq_ans[[1]], cortab=sq97)
O99seq <- A2O(ASVtable=seq_ans[[1]], cortab=sq99)

###### return genus, family table#####
ans <- taxaratio(seq_ans[[1]],seq_ans[[3]])
ans2 <- list()
for(i in 1:length(ans)){
  ans2[[i]] <- ans[[i]]
}
ans2[[(length(ans)+1)]] <- seq_ans[[1]]
ans2[[length(ans2)+1]] <- O99seq
ans2[[length(ans2)+1]] <- O97seq

#-- set save dir
saveRDS(ans2,sprintf("%s/%s_seqtab_%s.rds",fd,"List",as.character(number)))
saveRDS(seq_ans[[1]], sprintf("%s/%s",fd,seq_ans[[2]])); saveRDS(seq_ans[[3]], sprintf("%s/%s",fd,seq_ans[[4]]))
