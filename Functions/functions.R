library(cluster)
library(cowplot)
library(doParallel)
library(plyr)
library(dplyr)
library(ecoregime)
library(foreach)
library(ggforce)
library(ggnewscale)
library(ggplot2)
library(ggplate)
library(ggsci)
library(ggtext)
library(ggh4x)
library(gtools)
library(grid)
library(gridExtra)
library(igraph)
library(lemon)
library(lme4)
library(lmerTest)
library(lmodel2)
library(markdown)
library(MASS)
library(MCMCpack)
library(mgcv)
library(multimode)
library(openxlsx)
library(parallel)
library(philentropy)
library(plotly)
library(plyr)
library(ppcor)
library(RColorBrewer)
library(RasperGade)
library(RasperGade16S)
library(reshape2)
library(rELA)
library(scales)
library(seqinr)
library(stringdist)
library(stringr)
library(smacof)
library(tibble)
library(tidyr)
library(tidygraph)
library(tidyverse)
library(vegan)

taxaratio <- function(seqtab,taxa){
  Rows <- nrow(seqtab)
  f <- function(seqtab,taxa,taxid){ #- taxa=taxtable; taxid="Order"
    taxa_n <- cbind(ASV=rownames(taxa),ta=taxa[,(colnames(taxa)==taxid)])
    ls <- unique(taxa_n[,2]); ans <- matrix(0,nrow=nrow(seqtab),ncol=length(ls))
    for(i in 1:length(ls)){ #- i=1
      key <- taxa_n[(taxa_n[,2] == ls[i]),1]
      if(Rows>1){
        if(length(key)>1){
          ans[,i] <- rowSums(seqtab[,(colnames(seqtab) %in% key)])
        }else{
          ans[,i] <- seqtab[,(colnames(seqtab) %in% key)]
        }
      }else{
        if(length(key)>1){
          ans[,i] <- sum(seqtab[1,(colnames(seqtab) %in% key)])
        }else{
          ans[,i] <- seqtab[1,(colnames(seqtab) %in% key)]
        }
      }
    }
    rownames(ans) <- rownames(seqtab); colnames(ans) <- ls
    return(ans)
  }
  
  cla <- f(seqtab=seqtab,taxa=taxa,taxid="Class")
  ord <- f(seqtab=seqtab,taxa=taxa,taxid="Order")
  fam <- f(seqtab=seqtab,taxa=taxa,taxid="Family")
  gen <- f(seqtab=seqtab,taxa=taxa,taxid="Genus")
  
  return(list(cla,ord,fam,gen))
}

##############
palette_30 <- c(
  "#E6195B", "#3CB44B", "#FFE119", "#0082C8", "#F58231",
  "#911EB4", "#46F0F0", "#F032E6", "#D2F53C", "#FABEBE",
  "#008080", "#E6BEFF", "#AA6E28", "#FFFAC8", "#800000",
  "#AAFFC3", "#808000", "#FFD8B1", "#000080", "#AA8080",
  "#FFFFFF", "#0000AA", "#A9A9A9", "#B0E0E6", "#DC140C",
  "#ADFF2F", "#20B2AA", "#9370DB", "#FA8072", "#FF8C00"
)


# num is the number of taxa to be kept, and the rest will be merged into "Others"
sum_other <- function(mat,num){
  a <- colSums(mat)
  b <- sort(a, decreasing=TRUE)
  c <- (a>b[num])|(a==b[num])
  out <- matrix(0, nrow=nrow(mat), ncol=(sum(c)+1))
  count <- 1
  outname <- c(rep(NA,sum(c)),"Others")
  for(i in 1:ncol(mat)){
    if(c[i] == TRUE){
      out[,count] <- mat[,i]
      outname[count] <- colnames(mat)[i]
      count <- count+1
    }
    else{
      out[,num+1] <- out[,num+1]+mat[,i]
    }
  }
  rownames(out) <- rownames(mat)
  colnames(out) <- outname
  return(out)
}

rare_merge <- function(mat,csv,rarefuction){
  a <- rrarefy(mat, rarefuction)
  b <- cbind(rownames(a),as.data.frame(a))
  colnames(b) <- c("sample_name",colnames(b))
  c <- merge(b,csv,by="sample_name")
  
  return(c)
}

###############
lib <- 'RColorBrewer';library(package = lib, character.only=TRUE);packageVersion(lib) 
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col.vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

palettes <- function(x){ col.vector[1:length(unique(x))]}


###############
remove.outliers <- function(x, conf.level = 0.95)
{
  x <- x[!is.na(x)]
  del.val <- NULL
  
  while (TRUE) {
    n <- length(x)
    if (n < 3) {
      break
    }
    
    r <- range(x)
    t <- abs(r - mean(x)) / sd(x)
    q2 <- (n - 2) / ((n - 1) ^ 2 / t ^ 2 / n - 1)
    q2[q2 < 0] <- 0
    q <- sqrt(q2)
    p <- n * pt(q, n - 2, lower.tail = FALSE)
    
    if (t[1] < t[2]) {
      if (p[2] < 1 - conf.level) {
        del.val <- c(del.val, r[2])
        x <- x[x != r[2]]
        next
      }
    } else {
      if (p[1] < 1 - conf.level) {
        del.val <- c(del.val, r[1])
        x <- x[x != r[1]]
        next
      }
    }
    break
  }
  return(list(x = x, del.val = del.val))
}

forjac <- function(mat,thr=1,dataframe=FALSE){
  ans <- mat
  for(i in 1:nrow(ans)){
    for(j in 1:ncol(ans)){
      if(ans[i,j]>(thr-1)){
        ans[i,j] <- 1
      }else{
        ans[i,j] <- 0
      }
    }
  }
  if(dataframe==TRUE){
    ans <- as.data.frame(ans)
  }
  return(ans)
}

writefastafile=function(seqlist,namelist=seqlist,outfn="out"){
  #- Checking the arguments
  isbadarg=c(FALSE, "")
  if(!is.vector(namelist))isbadarg=c(TRUE, "argument ‘namelist’ is wrong type.")
  else if(!is.vector(seqlist))isbadarg=c(TRUE, "extra argument ‘seqlist’ is wrong type.")
  else if(!is.character(outfn))isbadarg=c(TRUE, "extra argument ‘outfn’ is wrong type.")
  else if(length(namelist)!=length(seqlist))isbadarg=c(TRUE, "Difference in vector length between seqlist and namelist")
  if(isbadarg[1]){
    warning(isbadarg[2])
    return()
  }
  library(seqinr)
  
  for(i in 1:length(namelist)){
    if(i == 1) mode = 'w'
    else mode = 'a'
    write.fasta(sequences = as.vector(seqlist[i]),names=as.vector(namelist[i]),file.out = paste(outfn, ".fasta",sep = ''), open=mode)
  }
}

Dictio <- function(vec, v1, v2){ #- vec=c(1,1,1,2,3,3,3,4); v1=c(1,2,3,4); v2=c("A","B","C","D")
  vec2 <- rep(NA, length(vec))
  for(i in 1:length(v1)){ #- i=1
    vec2[vec == v1[i]] <- v2[i]
  }
  return(vec2)
}
