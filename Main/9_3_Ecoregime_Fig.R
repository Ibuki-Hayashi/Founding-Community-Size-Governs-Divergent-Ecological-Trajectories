rm(list = ls(all.names = TRUE))

library(ggplot2); library(vegan);library(gridExtra);library(lemon);library(tidyr);library(markdown)
library(stringr);library(ggforce);library(cowplot);library(ggnewscale); library(reshape2);library(ggtext)
library(ecoregime)
source("Functions/functions.R")

fd <- "Ecoregime_Fig"

folder <- c(list.files("Output"), "hoge")
if(any(!(folder %in% fd))){
  dir.create(sprintf("Output/%s", fd))
}
fd <- sprintf("Output/%s", fd)

Daycol <- c("#1f78b4", "#33a02c", "#ff7f00", "#6a3d9a")
Inocol <- c("#4B0082", "#800080", "#C71585", "#E6CFE6")
SW <- c("#8B2500", "#004EbB")

# -- read data from previous dir
dfs_size <- read.csv("Output/Ecoregime_soil/figures/Minsegs_5/Metrices_OTU99_Size.csv")
dfw_size <- read.csv("Output/Ecoregime_water/figures/Minsegs_5/Metrices_OTU99_Size.csv")

dfs_link <- read.csv("Output/Ecoregime_soil/figures/Minsegs_5/Metrices_OTU99_Link.csv")
dfw_link <- read.csv("Output/Ecoregime_water/figures/Minsegs_5/Metrices_OTU99_Link.csv")

dfs_density <- read.csv("Output/Ecoregime_soil/figures/Minsegs_5/Metrices_OTU99_Density.csv")
dfw_density <- read.csv("Output/Ecoregime_water/figures/Minsegs_5/Metrices_OTU99_Density.csv")

dfs_depth <- read.csv("Output/Ecoregime_soil/figures/Minsegs_5/Metrices_OTU99_Depth.csv")
dfw_depth <- read.csv("Output/Ecoregime_water/figures/Minsegs_5/Metrices_OTU99_Depth.csv")

dfs_rank <- read.csv("Output/Ecoregime_soil/figures/Minsegs_5/Metrices_OTU99_Rank.csv")
dfw_rank <- read.csv("Output/Ecoregime_water/figures/Minsegs_5/Metrices_OTU99_Rank.csv")

FFF <- function(dfs, dfw){
  fd <- "Ecoregime_Fig"
  fd <- sprintf("Output/%s", fd)
  
  Daycol <- c("#1f78b4", "#33a02c", "#ff7f00", "#6a3d9a")
  Inocol <- c("#4B0082", "#800080", "#C71585", "#E6CFE6")
  SW <- c("#8B2500", "#004EbB")
  
  Density <- c(1, 0.1, 0.01, 0.001)
  
  df <- data.frame(Source=c(rep("Soil", 3), rep("Water", 3)), rbind(dfs, dfw))
  ldf <- gather(df, "Density", "Score", -Source, -X)
  ldf$Density <- factor(str_sub(ldf$Density, 9, -1), levels=Density)
  
  deve <- ldf[ldf$X == "dEv", ]
  ddis <- ldf[ldf$X == "dDis", ]
  dbd <- ldf[ldf$X == "dBD", ]
  
  peve <- ggplot(deve)+
    geom_line(aes(x=Density, y=Score, group=Source, linetype = Source), size=0.35)+
    geom_point(aes(x=Density, y=Score, fill=Density, shape=Source), size=1.8, color="grey10", stroke=0.4)+
    scale_fill_manual(values=Inocol)+theme_bw()+scale_shape_manual(values=c(21, 24))+
    xlab("Inoculum density")+ylab("Dynamic Evenness")+theme(text=element_text(size=7))
  pdis <- ggplot(ddis)+
    geom_line(aes(x=Density, y=Score, group=Source, linetype = Source), size=0.35)+
    geom_point(aes(x=Density, y=Score, fill=Density, shape=Source), size=1.8, color="grey10", stroke=0.4)+
    scale_fill_manual(values=Inocol)+theme_bw()+scale_shape_manual(values=c(21, 24))+
    xlab("Inoculum density")+ylab("Dynamic Dispersion")+theme(text=element_text(size=7), legend.position = "none")
  pbd <- ggplot(dbd)+
    geom_line(aes(x=Density, y=Score, group=Source, linetype = Source), size=0.35)+
    geom_point(aes(x=Density, y=Score, fill=Density, shape=Source), size=1.8, color="grey10", stroke=0.4)+
    scale_fill_manual(values=Inocol)+theme_bw()+scale_shape_manual(values=c(21, 24))+
    xlab("Inoculum density")+ylab("Dynamic Beta Diversity")+theme(text=element_text(size=7), legend.position = "none")
  
  ans <- plot_grid(pdis, pbd, peve, nrow=1, rel_widths=c(1, 1, 1.4))
  ans <- plot_grid(plot_grid(NA,NA,NA, nrow=1, labels=c("dDis", "dBD", "dEve"), label_size = 10), ans, nrow=2, rel_heights=c(0.1, 1))
  
  return(ans)
}

#-- Plot
pdf(sprintf("%s/Ecoregime_Metrics_Size.pdf", fd), width=8, height=2.2)
FFF(dfs=dfs_size, dfw=dfw_size)
dev.off()

pdf(sprintf("%s/Ecoregime_Metrics_Link.pdf", fd), width=8, height=2.2)
FFF(dfs=dfs_link, dfw=dfw_link)
dev.off()

pdf(sprintf("%s/Ecoregime_Metrics_Density.pdf", fd), width=8, height=2.2)
FFF(dfs=dfs_density, dfw=dfw_density)
dev.off()

pdf(sprintf("%s/Ecoregime_Metrics_Depth.pdf", fd), width=8, height=2.2)
FFF(dfs=dfs_depth, dfw=dfw_depth)
dev.off()

pdf(sprintf("%s/Ecoregime_Metrics_Rank.pdf", fd), width=8, height=2.2)
FFF(dfs=dfs_rank, dfw=dfw_rank)
dev.off()

