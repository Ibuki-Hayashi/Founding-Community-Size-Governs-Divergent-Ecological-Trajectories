source("Functions/functions.R")

#-- information of directory
fd <- "Output"; dd <- "Data"

#-- information of directory
exp <- read.csv(sprintf("%s/expdata_rep.csv", dd))
data <- readRDS(sprintf("%s/List_data_rep.rds", fd))
tax <- c("cla","ord","fam","gen","ASV","OTU99","OTU97")

InoJ <- readRDS(sprintf("%s/InoMulti_jac.rds", fd))[[6]]
Sil <- readRDS(sprintf("%s/Multinomialfir_Master_ACR_jac.rds", fd))[[6]]
Density <- c(1, 0.1, 0.01, 0.001)
day <- c(0, 2, 4, 6, 8)

Daycol <- c("#1f78b4", "#33a02c", "#ff7f00", "#6a3d9a")
Inocol <- c("#4B0082", "#800080", "#C71585", "#E6CFE6")

#-- 
InoJ$Day <- 0
dff <- as.data.frame(matrix(NA, ncol=ncol(Sil)-ncol(InoJ), nrow=nrow(InoJ)))
rownames(dff) <- rownames(InoJ)
colnames(dff) <- setdiff(colnames(Sil), colnames(InoJ))
InoJ2 <- cbind(InoJ, dff)

uSil <- bind_rows(InoJ2, Sil)

uSil$Day <- as.factor(uSil$Day)
uSil$Density <- factor(uSil$Density, levels=c(1, 0.1, 0.01, 0.001))
uSil$Statistic <- (uSil$Statistic-min(uSil$Statistic))/(max(uSil$Statistic)-min(uSil$Statistic))

uSilS <- uSil[uSil$Source == "Soil",]; uSilW <- uSil[uSil$Source == "Water",]

Dist_S <- ggplot(uSilS, aes(x=Day, y=Statistic))+
  xlab("Day")+ylab("Jaccard_multi")+theme_bw()+
  theme(text=element_text(size=7.2), legend.position="none", axis.title.y = element_markdown(size=7.2))+
  scale_color_manual(values=Inocol)+geom_line(aes(color=Density, group=Density))+scale_fill_manual(values=Inocol)+
  geom_point(aes(fill=Density), color="grey10", size=1.5, shape=21, stroke=0.3)

Dist_W <- ggplot(uSilW, aes(x=Day, y=Statistic))+
  xlab("Day")+ylab("Jaccard_multi")+theme_bw()+
  theme(text=element_text(size=7.2), legend.position="none", axis.title.y = element_markdown(size=7.2))+
  scale_color_manual(values=Inocol)+geom_line(aes(color=Density, group=Density))+scale_fill_manual(values=Inocol)+
  geom_point(aes(fill=Density), color="grey10", size=1.5, shape=24, stroke=0.3)

plpl <- plot_grid(Dist_S, Dist_W, ncol=2, nrow=1)

if(length(list.dirs(sprintf("%s/FigS18", fd), recursive=F))==0){
  dir.create(sprintf("%s/FigS18", fd))
}

pdf(sprintf("%s/FigS18/CommMultiDay0_8_Jaccard.pdf", fd), width=3.7, height=2)
plpl
dev.off()
