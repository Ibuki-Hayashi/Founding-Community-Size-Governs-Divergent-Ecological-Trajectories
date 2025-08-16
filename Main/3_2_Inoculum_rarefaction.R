rm(list = ls(all.names = TRUE))

library(ggplot2); library(vegan);library(gridExtra)
library(cowplot); library(markdown); library(ggtext)
library(tidyr); library(stringr);library(ggforce)
library(ggstar);library(ggnewscale);library(phyloseq)

source("Functions/functions.R")

#-- information of directory
fd <- "Inoculum"

folder <- c(list.files("Output"), "hoge")
if(any(!(folder %in% fd))){
  dir.create(sprintf("Output/%s", fd))
}
fd <- sprintf("Output/%s", fd)

dd1 <- "Data"; dd2 <- "Data2"

exp <- read.csv(sprintf("%s/expdata_rep.csv", dd1))
S_data <- readRDS(sprintf("%s/Inoculum_data_soil.rds", dd2))
W_data <- readRDS(sprintf("%s/Inoculum_data_water.rds", dd2))

sASV <- S_data[[1]][[5]]
wASV <- W_data[[1]][[5]]

# -- output
pdf(sprintf("%s/Inoculum_Rarefy_ASV.pdf", fd), width=8, height=4)
par(mfrow = c(1, 2))
rarecurve(x=t(sASV), step=50, label=T, ylim=c(0, 150), xlim=c(0, 33000),
          col="red", main="Soil-Inoculum", ylab="The number of ASVs", xlab="Reads")
rarecurve(t(wASV), step=50, label=T, ylim=c(0, 40), xlim=c(0, 4500),
          col="blue", main="Water-Inoculum", ylab="The number of ASVs", xlab="Reads")
dev.off()


