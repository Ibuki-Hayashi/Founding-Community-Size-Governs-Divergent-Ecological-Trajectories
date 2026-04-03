source("Functions/functions.R")

#-- information of directory
fd <- "Output"; dd <- "Data"

#-- information of directory
exp <- read.csv(sprintf("%s/expdata_rep.csv", dd))
data <- readRDS(sprintf("%s/List_data_rep.rds", fd))
Kuji <- readRDS(sprintf("%s/Inoculum_lottery_community.rds", fd))

tax <- c("cla","ord","fam","gen","ASV","OTU99","OTU97")
Density <- c(1, 0.1, 0.01, 0.001)
day <- c(2, 4, 6, 8)

Daycol <- c("#1f78b4", "#33a02c", "#ff7f00", "#6a3d9a")
Inocol <- c("#4B0082", "#800080", "#C71585", "#E6CFE6")

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
    kres <- pam(disdf, k=i, diss=TRUE)
    s <- silhouette(kres$cluster, dist(df))
    sil[i] <- mean(s[, 3])
  }
  ks <- c(1:6)[sil==max(sil)]; ans <- pam(disdf, k=ks, diss=TRUE)
  return(ans)
}

###################################
SilRess <- list()
qthr <- 0.05
Readthr <- 3000; fititr <- 10
for(i in 6){ #- i=6
  SilRes <- as.data.frame(matrix(NA, ncol=5, nrow=32))
  colnames(SilRes) <- c("Source", "Density", "Day", "Statistic", "pvalue")
  SilRes$Source <- c(rep("Soil", 16), rep("Water", 16))
  SilRes$Density <- rep(c(rep(1, 4), rep(0.1, 4), rep(0.01, 4), rep(0.001, 4)), 2)
  SilRes$Day <- rep(c(2, 4, 6, 8), 8)
  
  SilRes_HH <- SilRes_ACR <- SilRes

  for(j in 1:nrow(SilRes)){ #- j=1
    if(SilRes$Source[j] == "Soil"){
      ldf <- sTax[[i]][[which(Density == SilRes$Density[j])]][[which(day == SilRes$Day[j])]]
    }else{
      ldf <- wTax[[i]][[which(Density == SilRes$Density[j])]][[which(day == SilRes$Day[j])]]
    }
    ldf <- ldf[, c(7:ncol(ldf))]
    dis_ldf <- c(vegdist(ldf, method="bray"))
    
    Siltest <- modetest(data=dis_ldf, method="HH", B=5000)
    SilRes_HH$pvalue[j] <- Siltest$p.value
    SilRes_HH$Statistic[j] <- Siltest$statistic
    
    Siltest <- modetest(data=dis_ldf, method="ACR", B=5000)
    SilRes_ACR$pvalue[j] <- Siltest$p.value
    SilRes_ACR$Statistic[j] <- Siltest$statistic
  }
  
  Pgroup <- data.frame(ACR_p=SilRes_ACR$pvalue, HH_p=SilRes_HH$pvalue)
  ACRq <- (p.adjust(SilRes_ACR$pvalue, method="fdr") < 0.05)
  HHq <- (p.adjust(SilRes_HH$pvalue, method="fdr") < 0.05)
  QQ <- as.factor(ACRq + HHq)
  df <- data.frame(Pgroup, QQ=QQ)
  
  out1 <- ggplot(df, aes(x=ACR_p, y=HH_p, color=QQ))+
    geom_point(size=0.8)+theme_bw()+
    xlab("p-value (ACR method)")+ylab("p-value (HH method)")+
    geom_abline(slope=1, intercept=0, linetype="dashed", color="blue")+
    theme(text=element_text(size=10))+
    scale_color_manual(values=c("grey70", "orange", "#e31a1c"),
                       labels=c("Both Non-significant", "ACR or HH only", "Both significant"))
  
  df <- data.frame(ACR=SilRes_ACR$Statistic, HH=SilRes_HH$Statistic)
  
  out2 <- ggplot(df, aes(x=ACR, y=HH))+
    geom_point(size=0.3)+theme_bw()+
    xlab("Multimodality before standard (ACR method)")+ylab("Multimodality before standard (HH method)")+
    theme(text=element_text(size=10))
}

if(length(list.dirs(sprintf("%s/FigS6", fd), recursive=F))==0){
  dir.create(sprintf("%s/FigS6", fd))
}

pdf(sprintf("%s/FigS6/Comm_HHvsACR_CLR_pvalue.pdf", fd), width=4.5, height=3.5)
out1
dev.off()

pdf(sprintf("%s/FigS6/Comm_HHvsACR_CLR.pdf", fd), width=3.5, height=3.5)
out2
dev.off()
