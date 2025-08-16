rm(list = ls(all.names = TRUE))

library(ggplot2); library(vegan);library(gridExtra); library(cowplot);library(lemon); library(cluster)
library(tidyr); library(stringr);library(ggforce); library(ggtext);library(ggnewscale); library(ggh4x)
library(multimode); library(MCMCpack); library(MASS); library(reshape2); library(grid); library(openxlsx)

source("functions/functions.R")

#-- information of directory
fd <- "Quantification_Comm"

folder <- c(list.files("Output"), "hoge")
if(any(!(folder %in% fd))){
  dir.create(sprintf("Output/%s", fd))
}
fd <- sprintf("Output/%s", fd)

dd1 <- "Data"; dd2 <- "Data2"

exp <- read.csv(sprintf("%s/expdata_rep.csv", dd1))
data <- readRDS(sprintf("%s/List_data_rep.rds", dd2)); tax <- c("cla","ord","fam","gen","ASV","OTU99","OTU97")
Kuji <- readRDS(sprintf("%s/Inoculum_lottery_community.rds", dd2))
Sil <- readRDS(sprintf("%s/Multinomialfir_Master_ACR.rds", fd))

Sil2 <- list(); Sil2[[1]] <- Sil[[5]]; Sil2[[2]] <- Sil[[6]]; Sil2[[3]] <- Sil[[7]] 
names(Sil2) <- tax[5:7]
#-- Preservation as excel format
wb <- createWorkbook()
for (sheet_name in names(Sil2)) {
  addWorksheet(wb, sheet_name)  #- add a sheet
  writeData(wb, sheet_name, Sil2[[sheet_name]])  # wrting data
}
saveWorkbook(wb, sprintf("%s/Bray_Silres.xlsx", fd), overwrite = TRUE)
