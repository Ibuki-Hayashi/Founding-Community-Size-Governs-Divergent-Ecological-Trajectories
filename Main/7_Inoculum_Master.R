rm(list = ls(all.names = TRUE))

library(ggplot2); library(vegan);library(gridExtra); library(cowplot);library(lemon); library(cluster)
library(tidyr); library(stringr);library(ggforce); library(ggstar);library(ggnewscale); library(ggh4x)
library(multimode)

source("Functions/functions.R")

#-- information of directory
fd <- "Output"; dd <- "Data"
exp <- read.csv(sprintf("%s/expdata_rep.csv", dd))
data <- readRDS(sprintf("%s/List_data_rep.rds", fd)); tax <- c("cla","ord","fam","gen","ASV","OTU99","OTU97")
S_density <- c(1, 0.1, 0.01, 0.001)
W_density <- c(1, 0.1, 0.01, 0.001)
day <- c(2, 4, 6, 8)

Kuji <- readRDS(sprintf("%s/Inoculum_lottery.rds", fd))
sKuji <- Kuji[[1]]; wKuji <- Kuji[[2]]

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

#-- CV
CV_kuji <- function(ldf){ #- ldf <- sKuji[[4]][[5]]
  key <- unique(ldf$Taxonomy)
  ans <- numeric(length(key)); names(ans) <- key
  for(i in 1:length(key)){ #- i=10
    ans[i] <- sqrt(var(ldf$Copies[ldf$Taxonomy==key[i]]))/mean(ldf$Copies[ldf$Taxonomy==key[i]])
  }
  return(ans)
}

#-- CV caluculation
sCV <- list(); wCV <- list()
for(i in 5:7){ #- ASV, OTU99, OTU97
  sCV[[i]] <- list(); wCV[[i]] <- list()
  for(j in 1:4){ #- i=S_desnity, or W_density
    sCV[[i]][[j]] <- CV_kuji(sKuji[[j]][[i]])
    wCV[[i]][[j]] <- CV_kuji(wKuji[[j]][[i]])
  }
}

#-- Binomial distribution with normal error
log_likelihood <- function(sigma, observed_data, n, p){ #- n=3000; p=1/3; sigma=130; observed_data=df
  if (sigma <= 0) return(Inf)
  
  log_likelihood_value <- sum(log(sapply(observed_data, function(y) {
    sum(dbinom(0:n, size = n, prob = p) * dnorm(y, mean = 0:n, sd = sigma, log = FALSE))
  }, USE.NAMES = FALSE)))
  
  return(-log_likelihood_value)
}

#-- kmeans & Silhouette
ksil <- function(vec){ #- vec <- ldf
  vec2 <- vec+rnorm(length(vec), 0, 0.001)
  sil <- c(0)
  for(i in 2:6){ #- i=3
    kres <- kmeans(vec2, centers=i)
    s <- silhouette(kres$cluster, dist(vec2))
    sil[i] <- mean(s[, 3])  #- Average of Silhouette coefficient
  }
  k <- c(1:6)[sil==max(sil)]; ans <- kmeans(vec2, centers=k)
  return(ans)
}

#-- Robust scaling
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

#-- 1-1. Silverman test
SilRess <- list()
thr1 <- 0.005; thr2 <- 10 #-- thr1=relative abundance thr, thr2=occurence thr
qthr <- 0.05 #-- qthr = significance threshold of FDR
Readthr <- 3000; fititr <- 10 #-- Readthr = thr of Rarefaction、fititr = number of iteration of fitting
outl <- 5 #-- removed outliers

for(i in 5:7){ #- i=6
  SilRes <- as.data.frame(matrix(NA, ncol=10, nrow=0))
  colnames(SilRes) <- c("Source", "Density", "Day", "Taxa", "Count", "Mean", "pvalue",
                        "Statistic", "Sta_SD", "Sta_MAD")
  
  #-- Make df by Common taxa of each category
  for(j in 1:length(S_density)){ #- j=1
    for(k in 1:length(day)){ #- k=1
      ss <- sTax[[i]][[j]][[k]][,7:ncol(sTax[[i]][[j]][[k]])]; ws <- wTax[[i]][[j]][[k]][,7:ncol(wTax[[i]][[j]][[k]])]
      ss2 <- colnames(ss)[(colSums(ss)/sum(colSums(ss)) > thr1)&(sum(colSums(ss) > 0) > thr2)]
      ws2 <- colnames(ws)[(colSums(ws)/sum(colSums(ws)) > thr1)&(sum(colSums(ws) > 0) > thr2)]
      
      sCommon <- intersect(ss2, names(sCV[[i]][[j]]))
      wCommon <- intersect(ws2, names(wCV[[i]][[j]]))
      sMat <- data.frame(rep("Soil", length(sCommon)), rep(S_density[j], length(sCommon)),
                         rep(day[k], length(sCommon)), sCommon, matrix(NA, ncol=6, nrow=length(sCommon))); colnames(sMat) <- colnames(SilRes)
      wMat <- data.frame(rep("Water", length(wCommon)), rep(W_density[j], length(wCommon)),
                         rep(day[k], length(wCommon)), wCommon, matrix(NA, ncol=6, nrow=length(wCommon))); colnames(wMat) <- colnames(SilRes)
      
      SilRes <- rbind(SilRes, rbind(sMat, wMat))
    }
  }
  
  #-- Silverman test
  for(j in 1:nrow(SilRes)){
    #- Vector of used taxa
    if(SilRes$Source[j] == "Soil"){
      ldf <- sTax[[i]][[which(S_density == SilRes$Density[j])]][[which(day == SilRes$Day[j])]]
    }else{
      ldf <- wTax[[i]][[which(W_density == SilRes$Density[j])]][[which(day == SilRes$Day[j])]]
    }
    ldf <- as.vector(ldf[,colnames(ldf) == SilRes$Taxa[j]])
    SilRes$Count[j] <- sum(ldf>0)/length(ldf); SilRes$Mean[j] <- mean(ldf)
    
    ldf <- robust_scaling(ldf) #- scaling
    
    #- Silverman test
    Siltest <- modetest(data=ldf, method="ACR", B=3000)
    SilRes$pvalue[j] <- Siltest$p.value
    SilRes$Statistic[j] <- Siltest$statistic
    SilRes$Sta_SD[j] <- Siltest$statistic/sd(ldf)
    SilRes$Sta_MAD[j] <- Siltest$statistic/mad(ldf)
  }
  key <- sprintf("%s_%s_%s", SilRes$Source, SilRes$Density, SilRes$Day)
  
  for(sdd in 1:length(unique(key))){
    SilRes$qvalue[key==unique(key)[sdd]] <- p.adjust(SilRes$pvalue[key==unique(key)[sdd]], method="fdr")
  }
  
  #- Ino CV
  for(j in 1:nrow(SilRes)){ #- j=2
    if(SilRes$Source[j] == "Soil"){
      jj <- which(S_density == SilRes$Density[j])
      SilRes$Ino_CV[j] <- as.numeric(sCV[[i]][[jj]][SilRes$Taxa[j]])
    }else{
      jj <- which(W_density == SilRes$Density[j])
      SilRes$Ino_CV[j] <- as.numeric(wCV[[i]][[jj]][SilRes$Taxa[j]])
    }
  }
  
  #-- Fitting with Binomial distribution with normal error(for monomodal taxa only)
  AddSil <- data.frame(rep(NA, nrow(SilRes)), rep(NA, nrow(SilRes)), rep(NA, nrow(SilRes)))
  colnames(AddSil) <- c("Sigma1", "Loglikelihood1", "CV1")
  for(j in 1:nrow(SilRes)){ #- j=150
    #- Unimodal taxa
    if(SilRes$qvalue[j] > qthr){
      #- making using taxa
      if(SilRes$Source[j] == "Soil"){
        ldf <- sTax[[i]][[which(S_density == SilRes$Density[j])]][[which(day == SilRes$Day[j])]]
      }else{
        ldf <- wTax[[i]][[which(W_density == SilRes$Density[j])]][[which(day == SilRes$Day[j])]]
      }
      ldf <- as.vector(ldf[,colnames(ldf) == SilRes$Taxa[j]])
      ldf <- sort(ldf)[outl:(length(ldf)-outl)]
      
      binomp <- mean(ldf)/Readthr #- p in binomial dist(n = Readthr)
      #-- optimization
      itv <- c(); itll <- c()
      for(itr in 1:fititr){
        itv[itr] <- optimize(f=log_likelihood, interval=c(1, round(mean(ldf))+2), observed_data=ldf, n=Readthr, p=binomp)$minimum
        itll[itr] <- optimize(f=log_likelihood, interval=c(1, round(mean(ldf))+2), observed_data=ldf, n=Readthr, p=binomp)$objective
      }
      AddSil[j,1] <- mean(itv); AddSil[j,2] <- mean(itll); AddSil[j,3] <- mean(itv)/mean(ldf)
    }else{
      AddSil[j,1] <- NA; AddSil[j,2] <- NA; AddSil[j,3] <- NA
    }
  }
  SilRes <- cbind(SilRes, AddSil)
  
  #-- Fitting with Binomial distribution with normal error(for multimodal taxa only)
  AddSil <- data.frame(rep(NA, nrow(SilRes)), rep(NA, nrow(SilRes)), rep(NA, nrow(SilRes)))
  colnames(AddSil) <- c("Sigma2", "Loglikelihood2", "CV2")
  for(j in 1:nrow(SilRes)){ #- j=96
    #- unimodal taxa
    if(SilRes$qvalue[j] <= qthr){
      #- Vector of used taxa
      if(SilRes$Source[j] == "Soil"){
        ldf <- sTax[[i]][[which(S_density == SilRes$Density[j])]][[which(day == SilRes$Day[j])]]
      }else{
        ldf <- wTax[[i]][[which(W_density == SilRes$Density[j])]][[which(day == SilRes$Day[j])]]
      }
      ldf <- as.vector(ldf[,colnames(ldf) == SilRes$Taxa[j]])
      
      #- kmeans & Silhouette
      kres <- ksil(ldf)
      
      #- binomial fitting for each cluster
      clus <- unique(kres$cluster); itvM <- c(); itllM <- c(); itmM <- c()
      for(cl in 1:length(clus)){ #- cl=3
        ldf2 <- ldf[kres$cluster == clus[cl]]
        if((length(ldf2) > 5)&(var(ldf2)>0.1)){ #- Clusters including at least 6 taxa
          binomp <- mean(ldf2)/Readthr #- p in binomial dist(n = Readthr)
          #-- optimization
          itv <- c(); itll <- c()
          for(itr in 1:fititr){
            itv[itr] <- optimize(f=log_likelihood, interval=c(1, round(mean(ldf2))+2), observed_data=ldf2, n=Readthr, p=binomp)$minimum
            itll[itr] <- optimize(f=log_likelihood, interval=c(1, round(mean(ldf2))+2), observed_data=ldf2, n=Readthr, p=binomp)$objective
          }
          itvM <- c(itvM, mean(itv)); itllM <- c(itllM, mean(itll)); itmM <- c(itmM, mean(ldf2))
        }
      }

      AddSil[j,1] <- mean(itvM); AddSil[j,2] <- mean(itllM); AddSil[j,3] <- mean(itvM/itmM)
      
    }else{
      AddSil[j,1] <- NA; AddSil[j,2] <- NA; AddSil[j,3] <- NA
    }
  }
  
  SilRes <- cbind(SilRes, AddSil)
  
  SilRes[is.na(SilRes)] <- 0
  SilRes$sigma <- SilRes$Sigma1 + SilRes$Sigma2
  SilRes$LL <- SilRes$Loglikelihood1 + SilRes$Loglikelihood2
  SilRes$CV <- SilRes$CV1 + SilRes$CV2
  
  SilRess[[i]] <- SilRes
}

saveRDS(object=SilRess, file=sprintf("%s/Binomialfir_Master_HY.rds", fd))
