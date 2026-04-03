source("Functions/functions.R")

#-- information of directory
fd <- "Output"; dd <- "Data"

#-- information of directory
exp <- read.csv(sprintf("%s/expdata_rep.csv", dd))
data <- readRDS(sprintf("%s/List_data_rep.rds", fd))
tax <- c("cla","ord","fam","gen","ASV","OTU99","OTU97")

S_density <- c(1, 0.1, 0.01, 0.001)
W_density <- c(1, 0.1, 0.01, 0.001)
DEN <- c("x1", "x1/10", "x1/100", "x1/1000")
day <- c(2, 4, 6, 8)

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

#-- Loss species number
losssp <- function(mat){ #- mat=new_sKuji[[6]][[3]]
  all <- sum(colSums(mat)>0)
  ans <- c()
  for(i in 1:nrow(mat)){
    now <- sum(mat[i, ]>0)
    ans[i] <- all - now
  }
  return(ans)
}

#-- Panel A
da <- data.frame(Source=NA, Denity=NA, Losssp=NA)
for(i in 1:4){ #- i=1
  das <- losssp(new_sKuji[[6]][[i]]*3000/(rowSums(new_sKuji[[6]][[i]])[1]))
  daw <- losssp(new_wKuji[[6]][[i]]*3000/(rowSums(new_wKuji[[6]][[i]])[1]))
  da <- rbind(da, data.frame(Source=rep("Soil", length(das)), Denity=rep(DEN[i], length(das)), Losssp=das))
  da <- rbind(da, data.frame(Source=rep("Water", length(daw)), Denity=rep(DEN[i], length(daw)), Losssp=daw))
}
da <- da[-1, ]

pl_a <- ggplot(da)+
  geom_boxplot(aes(x=as.factor(Denity), y=Losssp, group=interaction(Source, as.factor(Denity))), 
               position=position_dodge(width=0.8), alpha=0.7, outlier.color = NA)+
  geom_jitter(aes(x=as.factor(Denity), y=Losssp, group=interaction(Source, as.factor(Denity))), 
                position=position_jitterdodge(jitter.width=0.2, dodge.width=0.8), size=0.3)+
  theme_bw()+
  xlab("Inoculum density")+
  ylab("Number of Extinct Species upon inoculation")+
  scale_x_discrete(labels=c("1"="1","0.1"="0.1","0.01"="0.01","0.001"="0.001"))+
  facet_wrap(~Source)+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=12),
        axis.title = element_text(size=14),
        axis.text = element_text(size=12))

if(length(list.dirs(sprintf("%s/FigS18", fd), recursive=F))==0){
  dir.create(sprintf("%s/FigS18", fd))
}

ggsave(sprintf("%s/FigS18/FigA_Inoculum_Extinction.pdf", fd), plot=pl_a, width=7, height=4)


#-- for Panel B
#-- Estimation the possibility of each taxon by Dirichlet distribution
Est_prob <- function(observed_counts, seed=1234){ #- ldf is count data (not relative abundance information)
  alpha_prior <- rep(1, length(observed_counts))
  alpha_posterior <- alpha_prior + observed_counts
  
  set.seed(seed)
  samples <- rdirichlet(1000, alpha_posterior)
  ans <- colMeans(samples); return(ans)
}

#-- kmeans & Silhouette
SilRess <- list()
qthr <- 0.05 #-- qthr
Readthr <- 3000; fititr <- 10 #-- Readthr
for(i in 6){ #- i=6
  SilRes <- as.data.frame(matrix(NA, ncol=5, nrow=8))
  colnames(SilRes) <- c("Source", "Density", "Day", "Statistic", "pvalue")
  SilRes$Source <- c(rep("Soil", 4), rep("Water", 4))
  SilRes$Density <- c(1, 0.1, 0.01, 0.001, 1, 0.1, 0.01, 0.001)

  #-- Silverman test
  for(j in 1:nrow(SilRes)){ #- j=1
    if(SilRes$Source[j] == "Soil"){
      ldf <- new_sKuji[[i]][[(1+log(SilRes$Density[j], 1/10))]]
    }else{
      ldf <- new_wKuji[[i]][[(1+log(SilRes$Density[j], 1/10))]]
    }
    dis_ldf <- c(vegdist(ldf, method="bray"))
    
    Siltest <- modetest(data=dis_ldf, method="ACR", B=50000)
    SilRes$pvalue[j] <- Siltest$p.value
    SilRes$Statistic[j] <- Siltest$statistic
  }
  SilRes$qvalue <- p.adjust(SilRes$pvalue, method="fdr")
  
  SilRess[[i]] <- SilRes
}

saveRDS(object=SilRess, file=sprintf("%s/InoMulti_bray.rds", fd))

SilRess <- list()
qthr <- 0.05 
Readthr <- 3000; fititr <- 10
for(i in 6){ #- i=6
  SilRes <- as.data.frame(matrix(NA, ncol=5, nrow=8))
  colnames(SilRes) <- c("Source", "Density", "Day", "Statistic", "pvalue")
  SilRes$Source <- c(rep("Soil", 4), rep("Water", 4))
  SilRes$Density <- c(1, 0.1, 0.01, 0.001, 1, 0.1, 0.01, 0.001)
  
  #-- Silverman test
  for(j in 1:nrow(SilRes)){ #- j=2
    if(SilRes$Source[j] == "Soil"){
      ldf <- new_sKuji[[i]][[(1+log(SilRes$Density[j], 1/10))]]
    }else{
      ldf <- new_wKuji[[i]][[(1+log(SilRes$Density[j], 1/10))]]
    }
    ldf <- forjac(ldf)
    dis_ldf <- c(vegdist(ldf, method="jaccard"))
    
    if(sum(dis_ldf)==0){
      SilRes$pvalue[j] <- 1
      SilRes$Statistic[j] <- 0
    }else{
      Siltest <- modetest(data=dis_ldf, method="ACR", B=50000)
      SilRes$pvalue[j] <- Siltest$p.value
      SilRes$Statistic[j] <- Siltest$statistic
    }
  }
  SilRes$qvalue <- p.adjust(SilRes$pvalue, method="fdr")
  
  SilRess[[i]] <- SilRes
}

saveRDS(object=SilRess, file=sprintf("%s/InoMulti_jac.rds", fd))
