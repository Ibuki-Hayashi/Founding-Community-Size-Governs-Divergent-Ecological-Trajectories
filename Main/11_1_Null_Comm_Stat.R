rm(list = ls(all.names = TRUE))

library(ggplot2); library(vegan);library(gridExtra); library(cowplot);library(lemon); library(cluster)
library(tidyr); library(stringr);library(ggforce); library(ggtext);library(ggnewscale); library(ggh4x)
library(multimode); library(MCMCpack); library(MASS); library(reshape2); library(grid); library(openxlsx)

source("Functions/functions.R")

#-- information of directory
fd <- "Output"; dd <- "Data"
Sil <- readRDS(sprintf("%s/Multinomialfir_Master_ACR.rds", fd))
tax <- c("cla","ord","fam","gen","ASV","OTU99","OTU97")

Sil2 <- list(); Sil2[[1]] <- Sil[[5]]; Sil2[[2]] <- Sil[[6]]; Sil2[[3]] <- Sil[[7]] 
names(Sil2) <- tax[5:7]
#-- Preservation as excel format
wb <- createWorkbook()
for (sheet_name in names(Sil2)) {
  addWorksheet(wb, sheet_name)  # シート追加
  writeData(wb, sheet_name, Sil2[[sheet_name]])  # データ書き込み
}
saveWorkbook(wb, sprintf("%s/Bray_Silres.xlsx", fd), overwrite = TRUE)
