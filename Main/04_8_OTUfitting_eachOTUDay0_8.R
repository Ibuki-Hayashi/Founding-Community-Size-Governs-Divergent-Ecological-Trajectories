source("Functions/functions.R")

#-- information of directory
fd <- "Output"; dd <- "Data"

#-- information of directory
exp <- read.csv(sprintf("%s/expdata_rep.csv", dd))
data <- readRDS(sprintf("%s/List_data_rep.rds", fd))
Sil <- readRDS(sprintf("%s/Binomialfir_Master_ACR.rds", fd))
inod <- readRDS(sprintf("%s/Inoculum_lottery_community.rds", fd))
tax <- c("cla","ord","fam","gen","ASV","OTU99","OTU97")

S_density <- c(1, 0.1, 0.01, 0.001)
W_density <- c(1, 0.1, 0.01, 0.001)
day <- c(2, 4, 6, 8)
Daycol <- c("#1f78b4", "#33a02c", "#ff7f00", "#6a3d9a")
Scol <- c("#8B2500", "#004EbB")
Inocol <- c("#4B0082", "#800080", "#C71585", "#E6CFE6")

#-- Function
robust_scaling <- function(x) {
  med <- median(x, na.rm = TRUE)
  iqr_val <- IQR(x, na.rm = TRUE)
  
  if (iqr_val == 0) {
    warning("IQR is 0, returning original data")
    return(x)
  }
  
  scaled_x <- (x - med) / iqr_val
  return(scaled_x)
}

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

#-- MA regression
plist <- list(); G_pl1 <- list(); G_pl2 <- list()
uden <- c("x1","x1/10","x1/100","x1/1000")

for(i in 6){ #- OTU99%
  #---------------
  ddf <- Sil[[i]]
  ddf$Density <- Dictio(ddf$Density, c(1,0.1,0.01,0.001), c("x1","x1/10","x1/100","x1/1000"))
  ddf$Density <- factor(ddf$Density, levels=c("x1","x1/10","x1/100","x1/1000"))
  ddf$Day <- factor(ddf$Day, levels=c(2,4,6,8))
  
  pl_m <- ggplot(ddf)+
    geom_boxplot(aes(x=Day, y=Statistic, group=Day), position=position_dodge(width=0.8), alpha=0.7, outlier.color = NA)+
    geom_jitter(aes(x=Day, y=Statistic, group=Day), position=position_jitterdodge(jitter.width=0.2, dodge.width=0.8), size=0.5)+
    theme_bw()+facet_grid(Source~as.factor(Density))
  
  pl_v <- ggplot(ddf)+
    geom_boxplot(aes(x=Day, y=CV, group=Day), position=position_dodge(width=0.8), alpha=0.7, outlier.color = NA)+
    geom_jitter(aes(x=Day, y=CV, group=Day), position=position_jitterdodge(jitter.width=0.2, dodge.width=0.8), size=0.5)+
    theme_bw()+facet_grid(Source~as.factor(Density))
  
  plist[[i]] <- plot_grid(pl_m, pl_v, nrow=2)
  
  #---------------
  plplp <- list()
  for(j in 1:length(S_density)){ #- j=1
    skeyotu <- list(); wkeyotu <- list()
    
    for(k in 1:length(day)){
      skeyotu[[k]] <- ddf$Taxa[(ddf$Source=="Soil")&(ddf$Density==uden[j])&(ddf$Day==day[k])]
      wkeyotu[[k]] <- ddf$Taxa[(ddf$Source=="Water")&(ddf$Density==uden[j])&(ddf$Day==day[k])]
    }
    skeyotu_all <- Reduce(intersect, skeyotu)
    wkeyotu_all <- Reduce(intersect, wkeyotu)
    
    sti <- paste("Soil", skeyotu_all, sep="_")
    wti <- paste("Water", wkeyotu_all, sep="_")
    
    spl <- list(); wpl <- list()
    for(k in 1:length(sti)){ #- k=1
      udff <- ddf[(ddf$Source=="Soil")&(ddf$Density==uden[j])&(ddf$Taxa==skeyotu_all[k]),]
      udff <- udff[, c("Day", "Statistic", "Taxa")]
      
      ino <- inod[[1]][[j]][[i]][[skeyotu_all[k]]]
      df0 <- data.frame(Day=0,
                        Statistic=modetest(data=robust_scaling(ino), method="ACR", B=5000)$statistic,
                        Taxa=skeyotu_all[k])
      udff <- rbind(df0, udff)
      udff$Sh <- as.factor(c(1,0,0,0,0))
      
      spl[[k]] <- ggplot(udff)+
        geom_point(aes(x=Day, y=Statistic, color=Sh, shape=Sh), size=1)+
        geom_line(aes(x=Day, y=Statistic, group=Taxa), size=0.4)+
        geom_hline(yintercept=udff[1,2], linetype="dashed", color="red", size=0.2)+
      
        scale_color_manual(values=c("1"="red", "0"="black"))+
        scale_shape_manual(values=c("1"=17, "0"=16))+
        theme_bw()+ylim(c(0, 0.35))+
        theme(
          axis.title = element_blank(),
          axis.text  = element_blank(),
          text = element_text(size=6.5),
          legend.position = "none"
        )+
        ggtitle(sti[k])
    }
    for(k in 1:length(wti)){
      udff <- ddf[(ddf$Source=="Water")&(ddf$Density==uden[j])&(ddf$Taxa==wkeyotu_all[k]),]
      udff <- udff[, c("Day", "Statistic", "Taxa")]
      
      ino <- inod[[2]][[j]][[i]][[wkeyotu_all[k]]]
      df0 <- data.frame(Day=0,
                        Statistic=modetest(data=robust_scaling(ino), method="ACR", B=5000)$statistic,
                        Taxa=wkeyotu_all[k])
      udff <- rbind(df0, udff)
      udff$Sh <- as.factor(c(1,0,0,0,0))
      
      wpl[[k]] <- ggplot(udff)+
        geom_point(aes(x=Day, y=Statistic, color=Sh, shape=Sh), size=1)+
        geom_line(aes(x=Day, y=Statistic, group=Taxa), size=0.4)+
        geom_hline(yintercept=udff[1,2], linetype="dashed", color="red", size=0.2)+
      
        scale_color_manual(values=c("1"="red", "0"="black"))+
        scale_shape_manual(values=c("1"=17, "0"=16))+
        theme_bw()+ylim(c(0, 0.35))+
        theme(
          axis.title = element_blank(),
          axis.text  = element_blank(),
          text = element_text(size=6.5),
          legend.position = "none"
        )+
        ggtitle(wti[k])
    }
    gws <- c(spl, wpl)
    if(length(gws) < 9){
      for(k in 1:(9-length(gws))){
        gws[[length(gws)+1]] <- ggplot() + theme_void()
      }
    }
    plplp[[j]] <- grid.arrange(grobs=gws, ncol=9)
  }
  G_pl1[[i]] <- grid.arrange(grobs=plplp, nrow=4)
  
  #-----
  plplp <- list()
  for(j in 1:length(S_density)){ #- j=1
    skeyotu <- list(); wkeyotu <- list()
    
    for(k in 1:length(day)){
      skeyotu[[k]] <- ddf$Taxa[(ddf$Source=="Soil")&(ddf$Density==uden[j])&(ddf$Day==day[k])]
      wkeyotu[[k]] <- ddf$Taxa[(ddf$Source=="Water")&(ddf$Density==uden[j])&(ddf$Day==day[k])]
    }
    skeyotu_all <- Reduce(intersect, skeyotu)
    wkeyotu_all <- Reduce(intersect, wkeyotu)
    
    sti <- paste("Soil", skeyotu_all, sep=" ")
    wti <- paste("Water", wkeyotu_all, sep=" ")
    
    spl <- list(); wpl <- list()
    for(k in 1:length(sti)){
      udff <- ddf[(ddf$Source=="Soil")&(ddf$Density==uden[j])&(ddf$Taxa==skeyotu_all[k]),]
      spl[[k]] <- ggplot(udff)+
        geom_point(aes(x=Day, y=CV), size=1)+
        geom_line(aes(x=Day, y=CV, group=Taxa), size=0.4)+
        theme(text = element_text(size=6))+
        theme_bw()+ylim(c(0, 1.5))+
        ggtitle(sti[k])
    }
    for(k in 1:length(wti)){
      udff <- ddf[(ddf$Source=="Water")&(ddf$Density==uden[j])&(ddf$Taxa==wkeyotu_all[k]),]
      wpl[[k]] <- ggplot(udff)+
        geom_point(aes(x=Day, y=CV), size=1)+
        geom_line(aes(x=Day, y=CV, group=Taxa), size=0.4)+
        theme(text = element_text(size=6))+
        theme_bw()+ylim(c(0, 1.5))+xlab("")+ylab("")+
        ggtitle(wti[k])
    }
    gws <- c(spl, wpl)
    if(length(gws) < 9){
      for(k in 1:(9-length(gws))){
        gws[[length(gws)+1]] <- ggplot() + theme_void()
      }
    }
    plplp[[j]] <- grid.arrange(grobs=gws, ncol=9)
  }
  G_pl2[[i]] <- grid.arrange(grobs=plplp, nrow=4)
}

if(length(list.dirs(sprintf("%s/FigS14", fd), recursive=F))==0){
  dir.create(sprintf("%s/FigS14", fd))
}

pdf(sprintf("%s/FigS14/EachOTU_Multimodality.pdf", fd), width=7.2, height=4)
plot(G_pl1[[6]])
dev.off()
