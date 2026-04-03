source("Functions/functions.R")

#-- information of directory
fd <- "Output"; dd <- "Data"

exp <- read.csv(sprintf("%s/expdata_rep.csv", dd))
S_data <- readRDS(sprintf("%s/Inoculum_data_soil.rds", fd))
W_data <- readRDS(sprintf("%s/Inoculum_data_water.rds", fd))

#- The list contains three objects (raw, copy, raspergade copy)
#- Copies = copies in 10uL

hoge <- S_data[[2]][[5]]
hoge <- hoge[hoge>1]

Kuji <- t(rmultinom(n=500, size=sum(hoge), prob=hoge/sum(hoge)))
k <- gather(as.data.frame(Kuji), key=Taxonomy, value=Copies)

sKuji <- list(); wKuji <- list()
sKuji_com <- list(); wKuji_com <- list()
for(j in 1:4){ #- j=Inoculum_density #- j=3
  sKuji[[j]] <- list(); wKuji[[j]] <- list()
  sKuji_com[[j]] <- list(); wKuji_com[[j]] <- list()
  for(i in 5:7){#- c(5:7)=c("ASV", "OTU99", "OTU97") #- i=5
    sdf <- S_data[[2]][[i]]/(10^(j-1))
    sdf <- sdf[sdf/sum(sdf) > 0.001]
    sKuji_com[[j]][[i]] <- as.data.frame(t(rmultinom(n=500, size=sum(sdf), prob=sdf/sum(sdf))))
    sKuji[[j]][[i]] <- gather(sKuji_com[[j]][[i]], key=Taxonomy, value=Copies)

    wdf <- W_data[[2]][[i]]/(10^(j-1))
    wdf <- wdf[wdf/sum(wdf) > 0.001]
    wKuji_com[[j]][[i]] <- as.data.frame(t(rmultinom(n=500, size=sum(wdf), prob=wdf/sum(wdf))))
    wKuji[[j]][[i]] <- gather(wKuji_com[[j]][[i]], key=Taxonomy, value=Copies)
  }
}

saveRDS(list(sKuji_com, wKuji_com), file=sprintf("%s/Inoculum_lottery_community.rds", fd))
saveRDS(list(sKuji, wKuji), file=sprintf("%s/Inoculum_lottery.rds", fd))

SAns_plot <- list(); WAns_plot <- list()
for(j in 1:4){
  SAns_plot[[j]] <- list(); WAns_plot[[j]] <- list()
  for(i in 5:7){
    SAns_plot[[j]][[i]] <- ggplot(sKuji[[j]][[i]])+
      geom_boxplot(aes(x=Taxonomy, y=Copies))+theme_bw()+
      theme(axis.text.x = element_text(angle=90))
    WAns_plot[[j]][[i]] <- ggplot(wKuji[[j]][[i]])+
      geom_boxplot(aes(x=Taxonomy, y=Copies))+theme_bw()+
      theme(axis.text.x = element_text(angle=90))
  }
}

if(length(list.dirs(sprintf("%s/FigS2_3_4", fd), recursive=F))==0){
  dir.create(sprintf("%s/FigS2_3_4", fd))
}

pdf(sprintf("%s/FigS2_3_4/Inoculum_lottery_OTU99_soil.pdf", fd),width=9,height=11)
plot(plot_grid(SAns_plot[[1]][[6]], SAns_plot[[2]][[6]], SAns_plot[[3]][[6]], SAns_plot[[4]][[6]], nrow=4))
dev.off()

pdf(sprintf("%s/FigS2_3_4/Inoculum_lottery_OTU99_water.pdf", fd),width=9,height=11)
plot(plot_grid(WAns_plot[[1]][[6]], WAns_plot[[2]][[6]], WAns_plot[[3]][[6]], WAns_plot[[4]][[6]], nrow=4))
dev.off()

##########################################
ssd <- data.frame(Dilution=NA, Taxonomy=NA, CV=NA, RA=NA)
wsd <- data.frame(Dilution=NA, Taxonomy=NA, CV=NA, RA=NA)
for(j in 1:4){ #- j=1
  tmp <- sKuji_com[[j]][[6]]
  CV <- apply(tmp, 2, function(x) sd(x)/mean(x))
  RA <- colMeans(tmp)/sum(colMeans(tmp))
  df <- data.frame(Dilution=rep(paste0("x", 10^(j-1)), length(CV)), Taxonomy=names(CV), CV=CV, RA=RA)
  ssd <- rbind(ssd, df)
  
  tmp <- wKuji_com[[j]][[6]]
  CV <- apply(tmp, 2, function(x) sd(x)/mean(x))
  RA <- colMeans(tmp)/sum(colMeans(tmp))
  df <- data.frame(Dilution=rep(paste0("x", 10^(j-1)), length(CV)), Taxonomy=names(CV), CV=CV, RA=RA)
  wsd <- rbind(wsd, df)
}
ssd <- ssd[-1, ]; wsd <- wsd[-1, ]
ssd2 <- ssd[ssd$RA>0.001,]; wsd2 <- wsd[wsd$RA>0.001,]
ssd2$RA <- log10(ssd2$RA); wsd2$RA <- log10(wsd2$RA)

p1 <- ggplot(ssd2, aes(x=RA, y=CV))+
  geom_point(size=0.5)+theme_bw()+
  labs(x="Log scale Relative abundance (%)", y="Coefficient of variation (CV)")+
  facet_wrap(~Dilution, ncol=4)+xlim(c(-3.1, -0.3))
p2 <- ggplot(wsd2, aes(x=RA, y=CV))+
  geom_point(size=0.5)+theme_bw()+
  labs(x="Log scale Relative abundance (%)", y="Coefficient of variation (CV)")+
  facet_wrap(~Dilution, ncol=4)+xlim(c(-3.1, -0.3))

pp <- plot_grid(p1, p2, nrow=2)

pdf(sprintf("%s/FigS2_3_4/CV_OTU99.pdf", fd), width=8.5, height=4.5)
plot(pp)
dev.off()
##########################################

ssd <- data.frame(Dilution=NA, Taxonomy=NA, CV=NA, RA=NA)
wsd <- data.frame(Dilution=NA, Taxonomy=NA, CV=NA, RA=NA)
for(j in 1:4){ #- j=1
  tmp <- sKuji_com[[j]][[6]]
  CV <- apply(tmp, 2, function(x) sd(x)/mean(x))
  RA <- colMeans(tmp)/sum(colMeans(tmp))
  df <- data.frame(Dilution=rep(paste0("x", 10^(j-1)), length(CV)), Taxonomy=names(CV), CV=CV, RA=RA)
  ssd <- rbind(ssd, df)
  
  tmp <- wKuji_com[[j]][[6]]
  CV <- apply(tmp, 2, function(x) sd(x)/mean(x))
  RA <- colMeans(tmp)/sum(colMeans(tmp))
  df <- data.frame(Dilution=rep(paste0("x", 10^(j-1)), length(CV)), Taxonomy=names(CV), CV=CV, RA=RA)
  wsd <- rbind(wsd, df)
}
ssd <- ssd[-1, ]; wsd <- wsd[-1, ]
ssd2 <- ssd[ssd$RA>0.001,]; wsd2 <- wsd[wsd$RA>0.001,]

p1 <- ggplot(ssd2, aes(x=RA, y=CV))+
  geom_point(size=0.5)+theme_bw()+
  labs(x="Log scale Relative abundance (%)", y="Coefficient of variation (CV)")+
  facet_wrap(~Dilution, ncol=4)#xlim(c(-3.1, -0.3))+
p2 <- ggplot(wsd2, aes(x=RA, y=CV))+
  geom_point(size=0.5)+theme_bw()+
  labs(x="Log scale Relative abundance (%)", y="Coefficient of variation (CV)")+
  facet_wrap(~Dilution, ncol=4)#xlim(c(-3.1, -0.3))+

pp <- plot_grid(p1, p2, nrow=2)

pdf(sprintf("%s/FigS2_3_4/CV_OTU99_notlog.pdf", fd), width=8.5, height=4.5)
plot(pp)
dev.off()
