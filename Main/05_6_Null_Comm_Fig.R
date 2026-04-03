source("Functions/functions.R")

#-- information of directory
fd <- "Output"; dd <- "Data"

#-- information of directory
exp <- read.csv(sprintf("%s/expdata_rep.csv", dd))
data <- readRDS(sprintf("%s/List_data_rep.rds", fd))
Kuji <- readRDS(sprintf("%s/Inoculum_lottery_community.rds", fd))
Sil <- readRDS(sprintf("%s/Multinomialfir_Master_ACR.rds", fd))

tax <- c("cla","ord","fam","gen","ASV","OTU99","OTU97")
Density <- c(1, 0.1, 0.01, 0.001)
day <- c(2, 4, 6, 8)

Daycol <- c("#1f78b4", "#33a02c", "#ff7f00", "#6a3d9a")
Inocol <- c("#4B0082", "#800080", "#C71585", "#E6CFE6")

#-- 
pl <- list(); plpl <- list()
for(i in 6){ #- i=6
  uSil <- Sil[[i]]
  uSil$Ino_CV <- log(uSil$Ino_CV, 10)
  uSil$Day <- as.factor(uSil$Day)
  uSil$Density <- factor(uSil$Density, levels=c(1, 0.1, 0.01, 0.001))
  uSil$sigma <- (uSil$sigma-min(uSil$sigma))/(max(uSil$sigma)-min(uSil$sigma))
  uSil$Statistic <- (uSil$Statistic-min(uSil$Statistic))/(max(uSil$Statistic)-min(uSil$Statistic))
  
  uSilS <- uSil[uSil$Source == "Soil",]; uSilW <- uSil[uSil$Source == "Water",]
  
  Cont_S <- ggplot(uSilS, aes(x=Day, y=sigma))+
    xlab("Day")+ylab("")+theme_bw()+
    ggtitle(sprintf("%s_Soil", tax[i]))+theme(text=element_text(size=7.2), legend.position="none", axis.title.y = element_markdown(size=7.2))+
    scale_color_manual(values=Inocol)+geom_line(aes(color=Density, group=Density))+scale_fill_manual(values=Inocol)+
    geom_point(aes(fill=Density), color="grey10", size=1.5, shape=21, stroke=0.3)
  Cont_W <- ggplot(uSilW, aes(x=Day, y=sigma))+
    xlab("Day")+ylab("")+theme_bw()+
    ggtitle(sprintf("%s_Water", tax[i]))+theme(text=element_text(size=7.2), legend.position="none", axis.title.y = element_markdown(size=7.2))+
    scale_color_manual(values=Inocol)+geom_line(aes(color=Density, group=Density))+scale_fill_manual(values=Inocol)+
    geom_point(aes(fill=Density), color="grey10", size=1.5, shape=24, stroke=0.3)
  
  Dist_S <- ggplot(uSilS, aes(x=Day, y=Statistic))+
    xlab("Day")+ylab("Multimodality in community-level")+theme_bw()+
    ggtitle(sprintf("%s_Soil", tax[i]))+theme(text=element_text(size=7.2), legend.position="none", axis.title.y = element_markdown(size=7.2))+
    scale_color_manual(values=Inocol)+geom_line(aes(color=Density, group=Density))+scale_fill_manual(values=Inocol)+
    geom_point(aes(fill=Density), color="grey10", size=1.5, shape=21, stroke=0.3)
  Dist_W <- ggplot(uSilW, aes(x=Day, y=Statistic))+
    xlab("Day")+ylab("Multimodality in community-level")+theme_bw()+
    ggtitle(sprintf("%s_Water", tax[i]))+theme(text=element_text(size=7.2), legend.position="none", axis.title.y = element_markdown(size=7.2))+
    scale_color_manual(values=Inocol)+geom_line(aes(color=Density, group=Density))+scale_fill_manual(values=Inocol)+
    geom_point(aes(fill=Density), color="grey10", size=1.5, shape=24, stroke=0.3)
  
  cl <- ggplot(uSil[uSil$Day == 8,])+
    geom_point(aes(x=Ino_CV, y=sigma, fill=Density), color="grey10", size=1.2, stroke=0.12)+theme_bw()+
    xlab("")+ylab("")+scale_fill_manual(values=Inocol)+
    ggtitle(sprintf("Day8"))+theme(text=element_text(size=7.2), axis.title.y = element_markdown(size=7.2))
  cl <- g_legend(cl)
  
  pl[[i]] <- plot_grid(plot_grid(Cont_S, Dist_S, Cont_W, Dist_W, ncol=2, nrow=2), cl, ncol=2, rel_widths=c(1, 0.1))
  plpl[[i]] <- plot_grid(plot_grid(Dist_S, Dist_W, Cont_S, Cont_W, ncol=4, nrow=1), cl, ncol=2, rel_widths=c(1, 0.1))
}

if(length(list.dirs(sprintf("%s/Fig3", fd), recursive=F))==0){
  dir.create(sprintf("%s/Fig3", fd))
}

pdf(sprintf("%s/Fig3/NullComm_plot.pdf", fd), width=8, height=7.2)
pl[[6]]
dev.off()

pdf(sprintf("%s/Fig3/NullComm_plot_forfig.pdf", fd), width=8, height=2)
plpl[[6]]
dev.off()
