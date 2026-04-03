source("Functions/functions.R")

#-- information of directory
fd <- "Output"; dd <- "Data"

#-- information of directory
exp <- read.csv(sprintf("%s/expdata_rep.csv", dd))
data <- readRDS(sprintf("%s/List_data_rep.rds", fd))
tax <- c("cla","ord","fam","gen","ASV","OTU99","OTU97")
Density <- c(1, 0.1, 0.01, 0.001)
day <- c(2, 4, 6, 8)

#-- data
sTax <- list(); wTax <- list()
for(i in 5:7){ #- i=1
  df <- data[[i]]; sdf <- df[df$Source=="soil",]; wdf <- df[df$Source=="water",]
  sTax[[i]] <- list(); wTax[[i]] <- list()
  for(j in 1:length(Density)){
    sTax[[i]][[j]] <- list(); wTax[[i]][[j]] <- list()
    for(k in 1:length(day)){
      sTax[[i]][[j]][[k]] <- sdf[(sdf$Inoculum_density == Density[j])&(sdf$Day == day[k]),]
      wTax[[i]][[j]][[k]] <- wdf[(wdf$Inoculum_density == Density[j])&(wdf$Day == day[k]),]
    }
  }
}

#-- 
Kuji <- readRDS(sprintf("%s/Inoculum_lottery_community.rds", fd))
sKuji <- Kuji[[1]]; wKuji <- Kuji[[2]]

new_sKuji <- list(); new_wKuji <- list()
for(i in 6){
  new_sKuji[[i]] <- list(); new_wKuji[[i]] <- list()
  for(j in 1:4){
    new_sKuji[[i]][[j]] <- sKuji[[j]][[i]][c(1:96),]
    new_wKuji[[i]][[j]] <- wKuji[[j]][[i]][c(1:96),]
  }
}
###################################

SilRes <- data.frame(Source=c(rep("Soil", 4), rep("Water", 4)),
                     Density=rep(Density[1:4], 2))
pls <- list()
for(i in 6){ #- i=6
  pls[[i]] <- list()
  for(jj in 1:8){ #- jj=4
    if(SilRes$Source[jj] == "Soil"){
      ldf <- sTax[[i]][[which(Density == SilRes$Density[jj])]][[which(day == 2)]]
    }else{
      ldf <- wTax[[i]][[which(Density == SilRes$Density[jj])]][[which(day == 2)]]
    }
    ldf <- ldf[, c(7:ncol(ldf))]
    ldf <- ldf[,(colSums(ldf)>0)]
    ldf <- forjac(as.matrix(ldf))
    d2_bc <- as.numeric(vegdist(ldf, method = "jaccard"))
    
    if(SilRes$Source[jj] == "Soil"){
      ldf <- new_sKuji[[i]][[which(Density == SilRes$Density[jj])]]
    }else{
      ldf <- new_wKuji[[i]][[which(Density == SilRes$Density[jj])]]
    }
    ldf <- ldf[,(colSums(ldf)>0)]
    ldf <- forjac(as.matrix(ldf))
    d0_bc <- as.numeric(vegdist(ldf, method = "jaccard"))
    
    pls[[i]][[jj]] <- local({
      Day0 <- d0_bc
      Day2 <- d2_bc
      ggplot()+
        geom_histogram(aes(x=Day2), bins=60, fill="green", alpha=0.7)+
        geom_histogram(aes(x=Day0), bins=60, fill="purple", alpha=0.7)+
        theme_bw()+theme(text=element_text(size=7),axis.title.y=element_markdown(),
                         legend.position="none")+xlab("Jaccard")+ylab("Frequency")
    })
  }
}

#-- Save
ans_soil <- plot_grid(pls[[6]][[1]], pls[[6]][[2]],
                      pls[[6]][[3]], pls[[6]][[4]], ncol=4, nrow=1)

ans_water <- plot_grid(pls[[6]][[5]], pls[[6]][[6]],
                       pls[[6]][[7]], pls[[6]][[8]], ncol=4, nrow=1)

if(length(list.dirs(sprintf("%s/FigS18", fd), recursive=F))==0){
  dir.create(sprintf("%s/FigS18", fd))
}

pdf(sprintf("%s/FigS18/Histo_Day02_soil.pdf", fd), width=8, height=1.3)
ans_soil
dev.off()

pdf(sprintf("%s/FigS18/Histo_Day02_water.pdf", fd), width=8, height=1.3)
ans_water
dev.off()
