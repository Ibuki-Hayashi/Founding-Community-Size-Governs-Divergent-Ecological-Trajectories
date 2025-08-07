if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("phyloseq")
BiocManager::install("treeio")
BiocManager::install("microbiome")

install.packages("ape")
install.packages("castor")
install.packages("vegan")
install.packages("seqinr")

if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")
devtools::install_github(repo = "wu-lab-uva/RasperGade")
devtools::install_github(repo = "wu-lab-uva/RasperGade16S", force=T)

library(vegan); library(stringr); library(RasperGade); library(RasperGade16S)

#-- information of directory
ans <- predict_16SGCN_from_sequences(seqs="table/merge_ASV_seq.fasta")

saveRDS(ans, file="HS3_eachASV.rds")
