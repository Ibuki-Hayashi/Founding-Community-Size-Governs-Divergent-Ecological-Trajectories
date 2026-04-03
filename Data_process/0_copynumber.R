library(vegan); library(stringr); library(RasperGade); library(RasperGade16S)

#-- information of directory
ans <- predict_16SGCN_from_sequences(seqs="Data/merge_ASV_seq.fasta")

saveRDS(ans, file="Output/HS3_eachASV.rds")
