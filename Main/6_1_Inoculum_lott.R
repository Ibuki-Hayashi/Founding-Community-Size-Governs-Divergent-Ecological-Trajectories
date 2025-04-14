rm(list = ls(all.names = TRUE))

library(ggplot2); library(vegan);library(gridExtra); library(cowplot);library(lemon); library(cluster)
library(tidyr); library(stringr);library(ggforce); library(ggstar);library(ggnewscale); library(ggh4x)

source("Functions/functions.R")

#-- information of directory
fd <- "Output"; dd <- "Data"

exp <- read.csv(sprintf("%s/expdata_rep.csv", dd))
S_data <- readRDS(sprintf("%s/Inoculum_data_soil.rds", fd))
W_data <- readRDS(sprintf("%s/Inoculum_data_water.rds", fd))

#- listが3つあるのは、それぞれ、生リード数、コピー数変換、RasperGade16S変換の結果
#- コピー数、RasperGade16Sはそれぞれ、10ul中に含まれる量

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
    sdf <- sdf[sdf/sum(sdf) > 0.000001]
    sKuji_com[[j]][[i]] <- as.data.frame(t(rmultinom(n=500, size=sum(sdf), prob=sdf/sum(sdf))))
    sKuji[[j]][[i]] <- gather(sKuji_com[[j]][[i]], key=Taxonomy, value=Copies)

    wdf <- W_data[[2]][[i]]/(10^(j-1))
    wdf <- wdf[wdf/sum(wdf) > 0.000001]
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

pdf(sprintf("%s/Inoculum_lottery_ASV_soil.pdf", fd),width=9,height=11)
plot(plot_grid(SAns_plot[[1]][[5]], SAns_plot[[2]][[5]], SAns_plot[[3]][[5]], SAns_plot[[4]][[5]], nrow=4))
dev.off()

pdf(sprintf("%s/Inoculum_lottery_OTU99_soil.pdf", fd),width=9,height=11)
plot(plot_grid(SAns_plot[[1]][[6]], SAns_plot[[2]][[6]], SAns_plot[[3]][[6]], SAns_plot[[4]][[6]], nrow=4))
dev.off()

pdf(sprintf("%s/Inoculum_lottery_OTU97_soil.pdf", fd),width=9,height=11)
plot(plot_grid(SAns_plot[[1]][[7]], SAns_plot[[2]][[7]], SAns_plot[[3]][[7]], SAns_plot[[4]][[7]], nrow=4))
dev.off()

pdf(sprintf("%s/Inoculum_lottery_ASV_water.pdf", fd),width=9,height=11)
plot(plot_grid(WAns_plot[[1]][[5]], WAns_plot[[2]][[5]], WAns_plot[[3]][[5]], WAns_plot[[4]][[5]], nrow=4))
dev.off()

pdf(sprintf("%s/Inoculum_lottery_OTU99_water.pdf", fd),width=9,height=11)
plot(plot_grid(WAns_plot[[1]][[6]], WAns_plot[[2]][[6]], WAns_plot[[3]][[6]], WAns_plot[[4]][[6]], nrow=4))
dev.off()

pdf(sprintf("%s/Inoculum_lottery_OTU97_water.pdf", fd),width=9,height=11)
plot(plot_grid(WAns_plot[[1]][[7]], WAns_plot[[2]][[7]], WAns_plot[[3]][[7]], WAns_plot[[4]][[7]], nrow=4))
dev.off()
