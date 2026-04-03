output.path <- "../Data"
## ======================================== ##

# -- If you want sum up sequence read,you select "sum".
type="sum"

referencepath="Ref/silva_nr99_v138_train_set_withSTD.fa" #-- please set your reference path here
vsearchpath="" #-- please set your vsearch path here
#-- 
##############################################

# -- Loading packages
library(dada2)
library(seqinr)

## -- Loading data table
files <- list.files()
tables=lapply(files, readRDS)

## -- Combine data table
stall <- mergeSequenceTables(tables=tables, repeats=type)
seqtab.nochim <- removeBimeraDenovo(stall, method="consensus", multithread=TRUE, verbose=TRUE)

taxa <- assignTaxonomy(seqtab.nochim, reference, multithread=TRUE)
taxa[is.na(taxa)] <- "Unidentified"

if(all(rownames(taxa)==colnames(seqtab.nochim))){
  
  seq.mat <- cbind(colnames(seqtab.nochim),sprintf('X_%s', formatC(1:ncol(seqtab.nochim), width = nchar(ncol(seqtab.nochim)), flag = "0")))
  write.fasta(as.list(seq.mat[,1]), seq.mat[,2], sprintf("%s/merge_ASV_seq.fasta", output.path) )

  colnames(seqtab.nochim) <- rownames(taxa) <- seq.mat[,2]
  ## ======================================== ##
  write.table(cbind(sample=rownames(seqtab.nochim), seqtab.nochim), sprintf("%s/merge_ASVtab.txt", output.path), sep="\t", quote=F, row.names=F) # CHANGE ME to where you want sequence table saved
  saveRDS(seqtab.nochim,  sprintf("%s/merge_ASVtab.rds", output.path))
  saveRDS(taxa,  sprintf("%s/merge_ASVtaxonomylist.rds", output.path))
  
  
}else{
  stop("rownames(taxa) and colnames(seqtab.nochim) are not match")
}

#####(OTU_99)##################################
minident <- 99
minident1 <- minident/100

system2(command = vsearchpath,
        args = c("--cluster_fast $input", sprintf("%s/merge_ASV_seq.fasta", output.path),
                 "--id", minident1,
                 "--mothur_shared_out", sprintf("%s/ASV_OTU_corestab_%s.txt", output.path, minident),
                 "--centroids", sprintf("%s/OTUseq_%s.fasta", output.path, minident),
                 "--msaout",  sprintf("%s/seqAlign_%s.txt", output.path, minident) ) )

otu <- read.table(sprintf("%s/ASV_OTU_corestab_%s.txt", output.path, minident),
                  header=TRUE, row.names=2)[,-c(1:2)]

## ========================================== ##
seqtab2 <- seqtab.nochim
otutab <- matrix(0, ncol=ncol(otu), nrow=nrow(seqtab2),
                 dimnames=list(rownames(seqtab2), colnames(otu)))

for(i in 1:ncol(otu)){
  
  if( sum(otu[,i])>1 ){
    memberSeq <- rownames(otu)[which(otu[,i]>0)]
    otutab[,i] <- rowSums(seqtab2[,which(colnames(seqtab2) %in% memberSeq) ])
  }else{
    centroidSeq <- colnames(otu)[i]
    otutab[,i] <- seqtab2[, which(colnames(seqtab2) == centroidSeq) ]
  }
  
}

saveRDS(otutab, sprintf("%s/seqOTUtab_%s.rds", output.path,minident))
write.csv(cbind(sample=rownames(otutab), otutab), sprintf('%s/seqOTUtab_%s.csv', output.path, minident), row.names=FALSE)
saveRDS(cbind(sample=rownames(otutab), otutab), sprintf('%s/seqOTUtab_%s.rds', output.path, minident))

#####(OTU_97)##################################
minident <- 97
minident1 <- minident/100

system2(command = vsearchpath,
        args = c("--cluster_fast $input", sprintf("%s/merge_ASV_seq.fasta", output.path),
                 "--id", minident1,
                 "--mothur_shared_out", sprintf("%s/ASV_OTU_corestab_%s.txt", output.path, minident),
                 "--centroids", sprintf("%s/OTUseq_%s.fasta", output.path, minident),
                 "--msaout",  sprintf("%s/seqAlign_%s.txt", output.path, minident) ) )

otu <- read.table(sprintf("%s/ASV_OTU_corestab_%s.txt", output.path, minident),
                  header=TRUE, row.names=2)[,-c(1:2)]

## ========================================== ##
seqtab2 <- seqtab.nochim
otutab <- matrix(0, ncol=ncol(otu), nrow=nrow(seqtab2),
                 dimnames=list(rownames(seqtab2), colnames(otu)))

for(i in 1:ncol(otu)){
  
  if( sum(otu[,i])>1 ){
    memberSeq <- rownames(otu)[which(otu[,i]>0)]
    otutab[,i] <- rowSums(seqtab2[,which(colnames(seqtab2) %in% memberSeq) ])
  }else{
    centroidSeq <- colnames(otu)[i]
    otutab[,i] <- seqtab2[, which(colnames(seqtab2) == centroidSeq) ]
  }
  
}

saveRDS(otutab, sprintf("%s/seqOTUtab_%s.rds", output.path,minident))
write.csv(cbind(sample=rownames(otutab), otutab), sprintf('%s/seqOTUtab_%s.csv', output.path, minident), row.names=FALSE)
saveRDS(cbind(sample=rownames(otutab), otutab), sprintf('%s/seqOTUtab_%s.rds', output.path, minident))
