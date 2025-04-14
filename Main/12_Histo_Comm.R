rm(list = ls(all.names = TRUE))

library(ggplot2); library(vegan);library(gridExtra); library(cowplot);library(lemon); library(cluster)
library(tidyr); library(stringr);library(ggforce); library(ggtext);library(ggnewscale); library(ggh4x)
library(multimode); library(MCMCpack); library(MASS); library(reshape2); library(grid)

source("Functions/functions.R")

#-- information of directory
fd <- "Output"; dd <- "Data"
exp <- read.csv(sprintf("%s/expdata_rep.csv", dd))
Kuji <- readRDS(sprintf("%s/Inoculum_lottery_community.rds", fd))
Sil <- readRDS(sprintf("%s/Multinomialfir_Master_ACR.rds", fd))
data <- readRDS(sprintf("%s/List_data_rep.rds", fd)); tax <- c("cla","ord","fam","gen","ASV","OTU99","OTU97")
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
    sil[i] <- mean(s[, 3])  # the average silhouette coefficient
  }
  ks <- c(1:6)[sil==max(sil)]; ans <- pam(disdf, k=ks, diss=T)
  return(ans)
}

###################################
SilRess <- list()
qthr <- 0.05 #-- qthr is the threshold of q-value for FDR
Readthr <- 3000; fititr <- 10 #-- Readthr is the thr of sequencing read、fititr is the number of iterations for the optimization
pls <- list()
for(i in 5:7){ #- i=6
  SilRes <- Sil[[i]]
  pls[[i]] <- list()
  for(jj in 1:nrow(SilRes)){ #- jj=2
    #- Unimodal dist
    if(SilRes$Source[jj] == "Soil"){
      ldf <- sTax[[i]][[which(Density == SilRes$Density[jj])]][[which(day == SilRes$Day[jj])]]
    }else{
      ldf <- wTax[[i]][[which(Density == SilRes$Density[jj])]][[which(day == SilRes$Day[jj])]]
    }
    ldf <- ldf[, c(7:ncol(ldf))]
    ldf <- ldf[,(colSums(ldf)>0)]
    if(SilRes$qvalue[jj] > qthr){
      size <- sum(ldf[1,])
      multi_probs <- Est_prob(colSums(ldf))
      n_samples <- nrow(ldf)
      n_categories <- ncol(ldf)
      
      itrs <- c(); itro <- c()
      for(itr in 1:100){
        fittt <- optimize(f=objective_function, interval=c(1, 10000), real_data=ldf)
        itrs[itr] <- fittt$minimum
        itro[itr] <- fittt$objective
      }
      sigma <- mean(itrs[(itro == min(itro))])
      
      sim_data <- t(rmultinom(n_samples, size = size, prob = multi_probs))
      var <- matrix(NA, ncol=n_categories, nrow=n_samples)
      for(ii in 1:ncol(var)){
        var[,ii] <- rnorm(n_samples, mean = 0, sd = sigma*multi_probs[ii])
      }
      sim_data <- sim_data + var
      sim_data[sim_data < 0] <- 0
      
      sim_bc <- as.numeric(vegdist(sim_data, method = "bray"))
      real_bc <- as.numeric(vegdist(ldf, method = "bray"))
      df <- data.frame(sim_bc, real_bc)
      
      pls[[i]][[jj]] <- ggplot(df)+
        geom_histogram(aes(x=real_bc), bins=60, fill="tomato", alpha=0.65)+
        geom_histogram(aes(x=sim_bc), bins=60, fill="lightblue", alpha=0.65)+
        theme_bw()+theme(text=element_text(size=7),axis.title.y=element_markdown(),
                         legend.position="none")+xlab("Bray-Curtis")+ylab("Frequency")+
        ggtitle(sprintf("%s_%s_%s", SilRes$Source[jj], SilRes$Day[jj], SilRes$Density[jj]))
    }else{
      kres <- ksil(ldf)
      clus <- unique(kres$cluster)
      sim <- matrix(NA, nrow=0, ncol=ncol(ldf))
      for(cl in 1:length(clus)){ #- cl=1
        ldf2 <- ldf[(kres$cluster == clus[cl]),]
        if(nrow(ldf2) > 5){ #- clusters with more than 5 samples
          size <- sum(ldf2[1,])
          multi_probs <- Est_prob(colSums(ldf2))
          n_samples <- nrow(ldf2)
          n_categories <- ncol(ldf2)
          
          itrs <- c(); itro <- c()
          for(itr in 1:100){
            fittt <- optimize(f=objective_function, interval=c(1, 10000), real_data=ldf2)
            itrs[itr] <- fittt$minimum
            itro[itr] <- fittt$objective
          }
          sigma <- mean(itrs[(itro == min(itro))])
          
          sim_data <- t(rmultinom(n_samples, size = size, prob = multi_probs))
          var <- matrix(NA, ncol=n_categories, nrow=n_samples)
          for(ii in 1:ncol(var)){
            var[,ii] <- rnorm(n_samples, mean = 0, sd = sigma*multi_probs[ii])
          }
          sim_data <- sim_data + var
          sim_data[sim_data < 0] <- 0
          
          sim <- rbind(sim, sim_data)
        }
      }
      if(nrow(sim) != nrow(ldf)){
        rep <- which(kres$id.med == min(kres$id.med[kres$id.med > 5]))
        if(length(rep) > 1){
          rep <- rep[1]
        }
        ldf2 <- ldf[(kres$cluster == clus[rep]),]
        size <- sum(ldf2[1,])
        multi_probs <- Est_prob(colSums(ldf2))
        n_samples <- nrow(ldf2)
        n_categories <- ncol(ldf2)
        
        itrs <- c(); itro <- c()
        for(itr in 1:100){
          fittt <- optimize(f=objective_function, interval=c(1, 10000), real_data=ldf2)
          itrs[itr] <- fittt$minimum
          itro[itr] <- fittt$objective
        }
        sigma <- mean(itrs[(itro == min(itro))])
        
        cor_samples <- nrow(ldf)-nrow(sim)
        sim_data <- t(rmultinom(cor_samples, size = size, prob = multi_probs))
        var <- matrix(NA, ncol=n_categories, nrow=cor_samples)
        for(ii in 1:ncol(var)){
          var[,ii] <- rnorm(cor_samples, mean = 0, sd = sigma*multi_probs[ii])
        }
        sim_data <- sim_data + var
        sim_data[sim_data < 0] <- 0
        
        sim <- rbind(sim, sim_data)
      }
      
      sim_bc <- as.numeric(vegdist(sim, method = "bray"))
      real_bc <- as.numeric(vegdist(ldf, method = "bray"))
      df <- data.frame(sim_bc, real_bc)
      
      pls[[i]][[jj]] <- ggplot(df)+
        geom_histogram(aes(x=real_bc), bins=60, fill="tomato", alpha=0.65)+
        geom_histogram(aes(x=sim_bc), bins=60, fill="lightblue", alpha=0.65)+
        theme_bw()+theme(text=element_text(size=7),axis.title.y=element_markdown(),
                         legend.position="none")+xlab("Bray-Curtis")+ylab("Frequency")+
        ggtitle(sprintf("%s_%s_%s", SilRes$Source[jj], SilRes$Day[jj], SilRes$Density[jj]))
    }
  }
}

#-- Save
ans_ASV <- grid.arrange(grobs=pls[[5]], ncol=5, nrow=7)
ans_99 <- grid.arrange(grobs=pls[[6]], ncol=5, nrow=7)
ans_97 <- grid.arrange(grobs=pls[[7]], ncol=5, nrow=7)

pdf(sprintf("%s/Histo_Bray_withModel.pdf", fd), width=8, height=11.5)
plot(ans_ASV)
plot(ans_99)
plot(ans_97)
dev.off()

pdf(sprintf("%s/Histo_Bray_withModel_OTU99.pdf", fd), width=8, height=11.5)
plot(ans_99)
dev.off()
