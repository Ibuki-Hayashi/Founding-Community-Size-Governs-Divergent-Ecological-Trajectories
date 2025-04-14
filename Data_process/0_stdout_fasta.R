rm(list = ls(all.names = TRUE))
library(seqinr)

#-- information of directory
fd <- "Output"

#-- read data from previous dir
source("Functions/functions.R")
fas <- read.alignment("Data/merge_ASV_seq.fasta", format="fasta")
taxtable <- readRDS("Data/merge_ASVtaxonomylist.rds")

##### fasta ###########
#- without stdDNA
name <- fas[[2]][taxtable[,3] != "c_"]
seqs <- fas[[3]][taxtable[,3] != "c_"]

writefastafile(seqs, namelist=name, outfn=sprintf("%s/ASV_outSTD", fd))

