rm(list = ls(all.names = TRUE))

library(ggplot2); library(vegan);library(gridExtra); library(cowplot);library(lemon); library(cluster)
library(tidyr); library(stringr);library(ggforce); library(ggstar);library(ggnewscale); library(ggh4x)
library(multimode); library(MCMCpack); library(doParallel)

source("functions/functions.R")

#-- information of directory
fd <- "Quantification_Comm"

folder <- c(list.files("Output"), "hoge")
if(any(!(folder %in% fd))){
  dir.create(sprintf("Output/%s", fd))
}
fd <- sprintf("Output/%s", fd)

dd1 <- "Data"; dd2 <- "Data2"

exp <- read.csv(sprintf("%s/expdata_rep.csv", dd1))
data <- readRDS(sprintf("%s/List_data_rep.rds", dd2)); tax <- c("cla","ord","fam","gen","ASV","OTU99","OTU97")

Kuji <- readRDS(sprintf("%s/Inoculum_lottery_community.rds", dd2))
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

#-- Optimizing function
objective_function <- function(sigma, real_data){ #- sigma=100
  size <- sum(real_data[1,])
  true_probs <- Est_prob(colSums(real_data))
  n_samples <- nrow(real_data)
  n_categories <- ncol(real_data)
  
  sim_data <- t(rmultinom(n_samples, size = size, prob = true_probs))
  var <- matrix(NA, ncol=ncol(real_data), nrow=n_samples)
  for(i in 1:ncol(var)){
    var[,i] <- rnorm(n_samples, mean = 0, sd = sigma*true_probs[i])
  }
  sim_data <- sim_data + var
  sim_data[sim_data < 0] <- 0
  sim_bc <- as.numeric(vegdist(sim_data, method = "bray"))
  real_bc <- as.numeric(vegdist(real_data, method = "bray"))
  
  ans <- mean((ecdf(sim_bc)(real_bc)-ecdf(real_bc)(sim_bc))^2)
  
  #ans <- sum((sim_bc - real_bc)^2)
  
  return(ans)
}

#-- Estimation the possibility of each taxon by Dirichlet distribution
Est_prob <- function(observed_counts, seed=1234){ #- ldf is count data (not relative abundance information)
  alpha_prior <- rep(1, length(observed_counts))
  alpha_posterior <- alpha_prior + observed_counts
  
  set.seed(seed)
  samples <- rdirichlet(1000, alpha_posterior)
  ans <- colMeans(samples); return(ans)
}

#-- kmeans & Silhouette
ksil <- function(df){ #- df <- ldf
  sil <- c(0)
  disdf <- vegdist(df, method="bray")
  for(i in 2:6){ #- i=3
    kres <- pam(disdf, k=i, diss=T)
    s <- silhouette(kres$cluster, dist(df))
    sil[i] <- mean(s[, 3])
  }
  ks <- c(1:6)[sil==max(sil)]; ans <- pam(disdf, k=ks, diss=T)
  return(ans)
}

###################################
SilRess <- list()
qthr <- 0.05
Readthr <- 3000; fititr <- 10
for(i in 5:7){ #- i=6
  SilRes <- as.data.frame(matrix(NA, ncol=5, nrow=32))
  colnames(SilRes) <- c("Source", "Density", "Day", "Statistic", "pvalue")
  SilRes$Source <- c(rep("Soil", 16), rep("Water", 16))
  SilRes$Density <- rep(c(rep(1, 4), rep(0.1, 4), rep(0.01, 4), rep(0.001, 4)), 2)
  SilRes$Day <- rep(c(2, 4, 6, 8), 8)
  
  #-- Silverman test
  for(j in 1:nrow(SilRes)){ #- j=1
    if(SilRes$Source[j] == "Soil"){
      ldf <- sTax[[i]][[which(Density == SilRes$Density[j])]][[which(day == SilRes$Day[j])]]
    }else{
      ldf <- wTax[[i]][[which(Density == SilRes$Density[j])]][[which(day == SilRes$Day[j])]]
    }
    ldf <- ldf[, c(7:ncol(ldf))]
    dis_ldf <- c(vegdist(ldf, method="bray"))
    
    #- Silverman test
    Siltest <- modetest(data=dis_ldf, method="ACR", B=50000)
    SilRes$pvalue[j] <- Siltest$p.value
    SilRes$Statistic[j] <- Siltest$statistic
  }
  SilRes$qvalue <- p.adjust(SilRes$pvalue, method="fdr")
  
  for(j in 1:nrow(SilRes)){ #- j=8
    if(SilRes$Source[j] == "Soil"){
      jj <- Kuji[[1]][[which(Density == SilRes$Density[j])]][[i]]
    }else{
      jj <- Kuji[[2]][[which(Density == SilRes$Density[j])]][[i]]
    }
    Var <- sum(diag(cov(t(jj))))/nrow(cov(t(jj)))
    CV <- sqrt(Var)/sum(jj[1,])
    SilRes$Ino_CV[j] <- CV
  }
  
  #-- Fitting with Binomial distribution with normal error(for monomodal taxa only)
  AddSil <- data.frame(rep(NA, nrow(SilRes)), rep(NA, nrow(SilRes)))
  colnames(AddSil) <- c("Sigma1", "ObFunc1")
  for(j in 1:nrow(SilRes)){ #- j=1
    #- Monomodal
    if(SilRes$qvalue[j] > qthr){
      if(SilRes$Source[j] == "Soil"){
        ldf <- sTax[[i]][[which(Density == SilRes$Density[j])]][[which(day == SilRes$Day[j])]]
      }else{
        ldf <- wTax[[i]][[which(Density == SilRes$Density[j])]][[which(day == SilRes$Day[j])]]
      }
      ldf <- ldf[, c(7:ncol(ldf))]
      ldf <- ldf[,(colSums(ldf)>0)]
      
      itrs <- c(); itro <- c()
      for(itr in 1:100){
        fittt <- optimize(f=objective_function, interval=c(1, 10000), real_data=ldf)
        itrs[itr] <- fittt$minimum
        itro[itr] <- fittt$objective
      }
      AddSil[j,1] <- mean(itrs[(itro == min(itro))]); AddSil[j,2] <- mean(itro[(itro == min(itro))])
      
    }else{
      AddSil[j,1] <- NA; AddSil[j,2] <- NA
    }
  }
  SilRes <- cbind(SilRes, AddSil)
  
  #-- Fitting with Binomial distribution with normal error(for multimodal taxa only)
  AddSil <- data.frame(rep(NA, nrow(SilRes)), rep(NA, nrow(SilRes)))
  colnames(AddSil) <- c("Sigma2", "ObFunc2")
  for(j in 1:nrow(SilRes)){ #- j=2
    #- Multi-modal
    if(SilRes$qvalue[j] <= qthr){
      if(SilRes$Source[j] == "Soil"){
        ldf <- sTax[[i]][[which(Density == SilRes$Density[j])]][[which(day == SilRes$Day[j])]]
      }else{
        ldf <- wTax[[i]][[which(Density == SilRes$Density[j])]][[which(day == SilRes$Day[j])]]
      }
      ldf <- ldf[, c(7:ncol(ldf))]
      ldf <- ldf[,(colSums(ldf)>0)]
      
      #- kmeans & Silhouette
      kres <- ksil(ldf)
      
      #- binomial fitting for each cluster
      clus <- unique(kres$cluster); sig <- c(); ff <- c()
      for(cl in 1:length(clus)){ #- cl=1
        ldf2 <- ldf[(kres$cluster == clus[cl]),]
        if(nrow(ldf2) > 5){
          itrs <- c(); itro <- c()
          for(itr in 1:100){
            fittt <- optimize(f=objective_function, interval=c(1, 10000), real_data=ldf2)
            itrs[itr] <- fittt$minimum
            itro[itr] <- fittt$objective
          }
          sig <- c(sig, mean(itrs[(itro == min(itro))])); ff <- c(ff, mean(itro[(itro == min(itro))]))
        }
      }
      AddSil[j,1] <- mean(sig); AddSil[j,2] <- mean(ff)
    }else{
      AddSil[j,1] <- NA; AddSil[j,2] <- NA
    }
  }
  
  SilRes <- cbind(SilRes, AddSil)
  
  SilRes[is.na(SilRes)] <- 0
  SilRes$sigma <- SilRes$Sigma1 + SilRes$Sigma2
  SilRes$ObFunc <- SilRes$ObFunc1 + SilRes$ObFunc2
  
  SilRess[[i]] <- SilRes
}

saveRDS(object=SilRess, file=sprintf("%s/Multinomialfir_Master_ACR.rds", fd))
