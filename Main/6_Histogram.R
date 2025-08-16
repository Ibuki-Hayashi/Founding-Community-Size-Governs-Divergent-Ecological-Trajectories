rm(list = ls(all.names = TRUE))

library(ggplot2); library(vegan);library(gridExtra);library(lemon);library(tidyr);library(markdown)
library(stringr);library(ggforce);library(cowplot);library(ggnewscale); library(reshape2);library(ggtext)
library(ecoregime); library(grid); library(cluster)
source("functions/functions.R")

fd <- "Histogram_OTU"
folder <- c(list.files("Output"), "hoge")
if(any(!(folder %in% fd))){
  dir.create(sprintf("Output/%s", fd))
}
fd <- sprintf("Output/%s", fd)

dd1 <- "Data"; dd2 <- "Data2"

#-- information of directory
exp <- read.csv(sprintf("%s/expdata_rep.csv", dd1))
data <- readRDS(sprintf("%s/List_data_rep.rds", dd2)); tax <- c("cla","ord","fam","gen","ASV","OTU99","OTU97")
Sil <- readRDS(sprintf("%s/Binomialfir_Master_HY.rds", dd2))

S_density <- c(1, 0.1, 0.01, 0.001)
W_density <- c(1, 0.1, 0.01, 0.001)
day <- c(2, 4, 6, 8)

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

#-- kmeans & Silhouette
ksil <- function(vec){ #- vec <- ldf
  vec2 <- vec+rnorm(length(vec), 0, 0.001)
  sil <- c(0)
  for(i in 2:6){ #- i=3
    kres <- kmeans(vec2, centers=i)
    s <- silhouette(kres$cluster, dist(vec2))
    sil[i] <- mean(s[, 3])  # Average of silhouette widths
  }
  k <- c(1:6)[sil==max(sil)]; ans <- kmeans(vec2, centers=k)
  return(ans)
}

log_likelihood <- function(sigma, observed_data, n, p){ #- n=3000; p=1/3; sigma=130; observed_data=df
  if (sigma <= 0) return(Inf)
  
  log_likelihood_value <- sum(log(sapply(observed_data, function(y) {
    sum(dbinom(0:n, size = n, prob = p) * dnorm(y, mean = 0:n, sd = sigma, log = FALSE))
  }, USE.NAMES = FALSE)))
  
  return(-log_likelihood_value)  # minus log-likelihood for minimization
}

#-- Histogram
qthr <- 0.05 #-- qthr is the threshold for q-value
Readthr <- 3000 #-- Readthr is the threshold for read count
histlist <- list()
for(i in 5:7){ #- i=5
  histlist[[i]] <- gList()
  SilRes <- Sil[[i]]
  for(j in 1:nrow(SilRes)){ #- j=4; j=260
    if(SilRes$Source[j] == "Soil"){
      ldf <- sTax[[i]][[which(S_density == SilRes$Density[j])]][[which(day == SilRes$Day[j])]]
    }else{
      ldf <- wTax[[i]][[which(W_density == SilRes$Density[j])]][[which(day == SilRes$Day[j])]]
    }
    ldf <- as.vector(ldf[,colnames(ldf) == SilRes$Taxa[j]])
    
    if(SilRes$qvalue[j] > qthr){
      x_vals <- seq(min(ldf), max(ldf)+300, length.out = 300) #- x-axis of the model
      var1 <- SilRes$Sigma1[j]
      n <- Readthr
      binomp <- mean(ldf)/n
      theory_density <- data.frame(
        x = x_vals,
        density = sapply(x_vals, function(x) sum(dbinom(0:n, size=n, prob=binomp) * dnorm(x - (0:n), mean=0, sd=var1)))
      )
      
      hist_info <- hist(c(ldf, min(ldf), max(ldf)+300), plot = FALSE, breaks = 60)
      
      theory_density$y_scaled <- theory_density$density * diff(hist_info$breaks[1:2]) * length(ldf)
      
      histlist[[i]][[j]] <- ggplot(data.frame(Sp=ldf), aes(x=Sp))+
        geom_histogram(bins = 60, fill = "lightblue", color = "black", alpha = 0.6, size=0.2)+
        theme_classic()+ggtitle(sprintf("%s_%s_%s_%s", SilRes$Source[j], SilRes$Day[j], SilRes$Density[j], SilRes$Taxa[j]))+
        theme(text=element_text(size=5.8),axis.title.y=element_markdown(),
              legend.position="none")+xlab("Relative abundance")+ylab("Frequency")+
        geom_line(data = theory_density, aes(x = x, y = y_scaled), color = "red", linewidth = 0.5)+
        scale_y_continuous(
          name = "Frequency",
          sec.axis = sec_axis(~ . / (diff(hist_info$breaks[1:2]) * length(ldf)),
                              name = "Density")
        )
    }else{
      #- kmeans & Silhouette
      kres <- ksil(ldf)
      
      #- binomial fitting for each cluster
      clus <- unique(kres$cluster); x_vals <- seq(min(ldf), max(ldf)+100, length.out = 300)
      density <- rep(0, length(x_vals))
      n <- Readthr
      for(cl in 1:length(clus)){ #- cl=2
        ldf2 <- ldf[kres$cluster == clus[cl]]
        if(length(ldf2) > 3){ #- Clusters which contains more than 3 data points
          binomp <- mean(ldf2)/Readthr #- p of binomial distribution
          #-- optimization(average of ten trials of optimization to avoid local minima)
          itv <- c()
          for(itr in 1:10){
            itv[itr] <- optimize(f=log_likelihood, interval=c(1, round(mean(ldf2))+2), observed_data=ldf2, n=Readthr, p=binomp)$minimum
          }
          itvM <- mean(itv)
          density1 <- sapply(x_vals, function(x) sum(dbinom(0:n, size=n, prob=binomp) * dnorm(x - (0:n), mean=0, sd=itvM)))
          density <- density + density1*(length(ldf2)/length(ldf))
        }
      }
      theory_density <- data.frame(
        x = x_vals,
        density = density
      )
      
      hist_info <- hist(c(ldf, min(ldf), max(ldf)+100), plot = FALSE, breaks = 60)
      
      theory_density$y_scaled <- theory_density$density * diff(hist_info$breaks[1:2]) * length(ldf)
      
      histlist[[i]][[j]] <- ggplot(data.frame(Sp=ldf), aes(x=Sp))+
        geom_histogram(bins = 60, fill = "lightblue", color = "black", alpha = 0.6, size=0.2)+
        theme_classic()+ggtitle(sprintf("%s_%s_%s_%s", SilRes$Source[j], SilRes$Day[j], SilRes$Density[j], SilRes$Taxa[j]))+
        theme(text=element_text(size=6),axis.title.y=element_markdown(),
              legend.position="none")+xlab("Relative abundance")+ylab("Frequency")+
        geom_line(data = theory_density, aes(x = x, y = y_scaled), color = "red", linewidth = 0.5)+
        scale_y_continuous(
          name = "Frequency",
          sec.axis = sec_axis(~ . / (diff(hist_info$breaks[1:2]) * length(ldf)),
                              name = "Density")
        )
    }
  }
}

cc <- gList()
for(i in 1:(240-217)){cc[[i]] <- nullGrob()}
ASV_p1 <- grid.arrange(grobs=histlist[[5]][c(1:40)], nrow=8, ncol=5)
ASV_p2 <- grid.arrange(grobs=histlist[[5]][c(41:80)], nrow=8, ncol=5)
ASV_p3 <- grid.arrange(grobs=histlist[[5]][c(81:120)], nrow=8, ncol=5)
ASV_p4 <- grid.arrange(grobs=histlist[[5]][c(121:160)], nrow=8, ncol=5)
ASV_p5 <- grid.arrange(grobs=histlist[[5]][c(161:200)], nrow=8, ncol=5)
ASV_p6 <- grid.arrange(grobs=c(c(histlist[[5]][c(201:217)]), cc), nrow=8, ncol=5)

pdf(sprintf("%s/ASV_histo_withModel.pdf", fd), width=8, height=11)
plot(ASV_p1)
plot(ASV_p2)
plot(ASV_p3)
plot(ASV_p4)
plot(ASV_p5)
plot(ASV_p6)
dev.off()

cc <- gList()
for(i in 1:(200-175)){cc[[i]] <- nullGrob()}
OTU99_p1 <- grid.arrange(grobs=histlist[[6]][c(1:40)], nrow=8, ncol=5)
OTU99_p2 <- grid.arrange(grobs=histlist[[6]][c(41:80)], nrow=8, ncol=5)
OTU99_p3 <- grid.arrange(grobs=histlist[[6]][c(81:120)], nrow=8, ncol=5)
OTU99_p4 <- grid.arrange(grobs=histlist[[6]][c(121:160)], nrow=8, ncol=5)
OTU99_p5 <- grid.arrange(grobs=c(c(histlist[[6]][c(161:175)]), cc), nrow=8, ncol=5)

pdf(sprintf("%s/OTU99_histo_withModel.pdf", fd), width=8, height=11)
plot(OTU99_p1)
plot(OTU99_p2)
plot(OTU99_p3)
plot(OTU99_p4)
plot(OTU99_p5)
dev.off()

#-- 99 devided
pdf(sprintf("%s/OTU99_s1.pdf", fd), width=8, height=11)
plot(OTU99_p1)
dev.off()
pdf(sprintf("%s/OTU99_s2.pdf", fd), width=8, height=11)
plot(OTU99_p2)
dev.off()
pdf(sprintf("%s/OTU99_s3.pdf", fd), width=8, height=11)
plot(OTU99_p3)
dev.off()
pdf(sprintf("%s/OTU99_s4.pdf", fd), width=8, height=11)
plot(OTU99_p4)
dev.off()
pdf(sprintf("%s/OTU99_s5.pdf", fd), width=8, height=11)
plot(OTU99_p5)
dev.off()


cc <- gList()
for(i in 1:(120-115)){cc[[i]] <- nullGrob()}
OTU97_p1 <- grid.arrange(grobs=histlist[[7]][c(1:40)], nrow=8, ncol=5)
OTU97_p2 <- grid.arrange(grobs=histlist[[7]][c(41:80)], nrow=8, ncol=5)
OTU97_p3 <- grid.arrange(grobs=c(c(histlist[[7]][c(81:115)]), cc), nrow=8, ncol=5)

pdf(sprintf("%s/OTU97_histo_withModel.pdf", fd), width=8, height=11)
plot(OTU97_p1)
plot(OTU97_p2)
plot(OTU97_p3)
dev.off()
