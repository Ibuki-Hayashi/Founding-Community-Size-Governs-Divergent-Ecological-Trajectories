rm(list = ls(all.names = TRUE))

library(ggplot2); library(vegan);library(gridExtra);library(lemon);library(tidyr);library(markdown)
library(stringr);library(ggforce);library(cowplot);library(ggnewscale); library(reshape2);library(ggtext)
library(ecoregime); library(grid)
source("Functions/functions.R")

fd <- "Output"; dd <- "Data"

#-- information of directory
exp <- read.csv(sprintf("%s/expdata_rep.csv", dd))
data <- readRDS(sprintf("%s/List_data_rep.rds", fd)); tax <- c("cla","ord","fam","gen","ASV","OTU99","OTU97")
Sil <- readRDS(sprintf("%s/Binomialfir_Master_HY.rds", fd))
S_density <- W_density <- c(1, 0.1, 0.01, 0.001)
day <- c(2, 4, 6, 8)

Daycol <- c("#1f78b4", "#33a02c", "#ff7f00", "#6a3d9a")
Inocol <- c("#4B0082", "#800080", "#C71585", "#E6CFE6")

#-- data
sTax <- list(); wTax <- list()
for(i in 5:7){ #- i=1
  df <- data[[i]]; sdf <- df[df$Source=="soil",]; wdf <- df[df$Source=="water",]
  sTax[[i]] <- list(); wTax[[i]] <- list()
  for(j in 1:4){
    sTax[[i]][[j]] <- list(); wTax[[i]][[j]] <- list()
    for(k in 1:length(day)){
      sTax[[i]][[j]][[k]] <- sdf[(sdf$Inoculum_density==S_density[j])&(sdf$Day == day[k]),]
      wTax[[i]][[j]][[k]] <- wdf[(wdf$Inoculum_density==W_density[j])&(wdf$Day == day[k]),]
    }
  }
}

#-- Silplot
corpl <- list(); Sanko <- list()
corpl1 <- list(); corpl2 <- list()
for(i in 6){ #- i=6
  SilRes <- Sil[[i]]
  
  uSil <- SilRes[(SilRes$LL != Inf),]
  uSil <- SilRes[(SilRes$CV != Inf),]
  uSil$Ino_CV <- log(uSil$Ino_CV, 10)
  uSil$unique <- sprintf("%s_%s", uSil$Source, uSil$Taxa)
  uSil$Day <- as.factor(uSil$Day)
  uSil$Density <- factor(uSil$Density, levels=c(1, 0.1, 0.01, 0.001))
  uSil$CV <- (uSil$CV-min(uSil$CV))/(max(uSil$CV)-min(uSil$CV))
  uSil$Statistic <- (uSil$Statistic-min(uSil$Statistic))/(max(uSil$Statistic)-min(uSil$Statistic))
  
  corpl1[[i]] <- gList(); corpl2[[i]] <- gList()
  itr <- 1
  for(j in 1:length(unique(uSil$unique))){ #- j=1
    df <- uSil[uSil$unique == unique(uSil$unique)[j],]
    key <- c()
    for(k in 1:4){
      key[k] <- length(df$Day[df$Day == day[k]])
    }
    if(sum(key > 2) >= 2){
      corpl1[[i]][[itr]] <- ggplot(df, aes(x=Ino_CV, y=CV, color=as.factor(Day)))+
        geom_point(size=1.3)+geom_line(size=0.25)+theme_bw()+xlab("Initial variation<br>upon inoculation")+ylab("Varation caused<br>after inoculation")+
        ggtitle(sprintf("%s", unique(uSil$unique)[j]))+theme(legend.position="none", text=element_text(size=6.5), axis.title=element_markdown())+
        scale_color_manual(values=Daycol)
      corpl2[[i]][[itr]] <- ggplot(df, aes(x=Ino_CV, y=Statistic, color=as.factor(Day)))+
        geom_point(size=1.3)+geom_line(size=0.25)+theme_bw()+xlab("Initial variation<br>upon inoculation")+ylab("Multimodality")+
        ggtitle(sprintf("%s", unique(uSil$unique)[j]))+theme(legend.position="none", text=element_text(size=6.5), axis.title=element_markdown())+
        scale_color_manual(values=Daycol)
      
      itr <- itr+1
    }
  }
  leg <- ggplot(data=data.frame(x=c(2,4,6,8), Day=as.factor(c(2,4,6,8))))+
    geom_point(aes(x=x, y=Day, color=Day))+theme_bw()+scale_color_manual(values=Daycol)+
    theme(text=element_text(size=6.5))
  leg <- g_legend(leg)
  
  n1 <- grid.arrange(grobs=corpl1[[i]], ncol=6, nrow=2)
  n1 <- plot_grid(n1, leg, ncol=2, rel_widths=c(1, 0.1))
  n2 <- grid.arrange(grobs=corpl2[[i]], ncol=6, nrow=2)
  n2 <- plot_grid(n2, leg, ncol=2, rel_widths=c(1, 0.1))
  corpl[[i]] <- plot_grid(n1, n2, ncol=1, rel_heights=c(1, 1))
  
  topCD <- ggplot(uSil)+
    geom_point(aes(x=CV, y=Statistic, color=Day), size=1.5)+theme_bw()+
    xlab("Varation caused after inoculation")+ylab("Multimodality")+
    ggtitle(sprintf("%s_St VS Mt", tax[i]))+theme(text=element_text(size=7.2), axis.title=element_markdown())+
    scale_color_manual(values=Daycol)
  
  cd2 <- ggplot(uSil[uSil$Day == 2,])+
    geom_point(aes(x=CV, y=Statistic, color=Density), size=0.8)+theme_bw()+
    xlab("Varation caused after inoculation")+ylab("Multimodality")+scale_color_manual(values=Inocol)+
    ggtitle(sprintf("%s_Day2", tax[i]))+theme(text=element_text(size=7.2), legend.position = "none", axis.title=element_markdown())
  cd4 <- ggplot(uSil[uSil$Day == 4,])+
    geom_point(aes(x=CV, y=Statistic, color=Density), size=0.8)+theme_bw()+
    xlab("Varation caused after inoculation")+ylab("Multimodality")+scale_color_manual(values=Inocol)+
    ggtitle(sprintf("%s_Day4", tax[i]))+theme(text=element_text(size=7.2), legend.position = "none", axis.title=element_markdown())
  cd6 <- ggplot(uSil[uSil$Day == 6,])+
    geom_point(aes(x=CV, y=Statistic, color=Density), size=0.8)+theme_bw()+
    xlab("Varation caused after inoculation")+ylab("Multimodality")+scale_color_manual(values=Inocol)+
    ggtitle(sprintf("%s_Day6", tax[i]))+theme(text=element_text(size=7.2), legend.position = "none", axis.title=element_markdown())
  cd8 <- ggplot(uSil[uSil$Day == 8,])+
    geom_point(aes(x=CV, y=Statistic, color=Density), size=0.8)+theme_bw()+
    xlab("Varation caused after inoculation")+ylab("Multimodality")+scale_color_manual(values=Inocol)+
    ggtitle(sprintf("%s_Day8", tax[i]))+theme(text=element_text(size=7.2), legend.position = "none", axis.title=element_markdown())
  
  leg <- ggplot(data=data.frame(x=c(2,4,6,8), Density=factor(c(1,0.1,0.01,0.001), levels=c(1,0.1,0.01,0.001))) )+
    geom_point(aes(x=x, y=Density, color=Density))+theme_bw()+scale_color_manual(values=Inocol)+
    theme(text=element_text(size=7.2))
  leg <- g_legend(leg)
  
  Sanko[[i]] <- plot_grid(topCD, plot_grid(plot_grid(cd2, cd4, cd6, cd8, nrow=2, ncol=2), leg, ncol=2, rel_widths = c(1,0.3)),
                          ncol=2, rel_widths = c(0.8,1))
}

Fig <- plot_grid(corpl[[6]], Sanko[[6]], ncol=1, rel_heights=c(1, 0.65),
                 labels=c("A", "B"), label_size=10)

#-- Save
pdf(sprintf("%s/InoCV_Cor_SFig.pdf", fd), width=8, height=9)
plot(Fig)
dev.off()
