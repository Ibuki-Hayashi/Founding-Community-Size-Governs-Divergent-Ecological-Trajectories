source("Functions/functions.R")

#-- information of directory
fd <- "Output"; dd <- "Data"

#-- information of directory
exp <- read.csv(sprintf("%s/expdata_rep.csv", dd))
data <- readRDS(sprintf("%s/List_data_rep.rds", fd))
tax <- c("cla","ord","fam","gen","ASV","OTU99","OTU97")

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

#-- CLR conversion
clr <- function(x, pseudocount = 1e-6) { #- x=a
  x <- x + pseudocount
  gm <- exp(rowMeans(log(x)))
  log(x / gm)
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
    sil[i] <- mean(s[, 3])  # the average silhouette width for the current number of clusters
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

#-- Multimode test by ACR and HH method fro each OTU (compare two methods)
SilRess <- list()
thr1 <- 0.005; thr2 <- 10 #-- thr1 is the threshold for relative abundance, thr2 is the threshold for the number of samples in which the taxon is present (both are used to determine the common taxa)
qthr <- 0.05 #-- qthr is the threshold for q-value to determine whether a taxon is monomodal or multimodal (monomodal if q-value > qthr, multimodal if q-value <= qthr)
Readthr <- 3000; fititr <- 10 #-- Readthr is the threshold of sequencing read, fititr is the number of iterations for optimization (the average of the results of fititr iterations is adopted as the final result)
outl <- 5 #-- the number of outliers to be removed from the data for fitting with binomial distribution with normal error (for monomodal taxa only)
otl <- list()

for(i in 6){ #- i=6
  SilRes <- as.data.frame(matrix(NA, ncol=10, nrow=0))
  colnames(SilRes) <- c("Source", "Density", "Day", "Taxa", "Count", "Mean", "pvalue",
                        "Statistic", "Sta_SD", "Sta_MAD")
  otl[[i]] <- as.data.frame(matrix(NA, ncol=6, nrow=0))
  colnames(otl[[i]]) <- c("Source", "Density", "Day", "Taxa", "Sample_name", "Wellpostition")
  
  #-- Make df by Common taxa of each category
  for(j in 1:length(S_density)){ #- j=1
    for(k in 1:length(day)){ #- k=1
      ss <- sTax[[i]][[j]][[k]][,7:ncol(sTax[[i]][[j]][[k]])]; ws <- wTax[[i]][[j]][[k]][,7:ncol(wTax[[i]][[j]][[k]])]
      ss2 <- colnames(ss)[(colSums(ss)/sum(colSums(ss)) > thr1)&(colSums(ss>0) > thr2)]
      ws2 <- colnames(ws)[(colSums(ws)/sum(colSums(ws)) > thr1)&(colSums(ws>0) > thr2)]
      
      sCommon <- intersect(ss2, names(sCV[[i]][[j]]))
      wCommon <- intersect(ws2, names(wCV[[i]][[j]]))
      sMat <- data.frame(rep("Soil", length(sCommon)), rep(S_density[j], length(sCommon)),
                         rep(day[k], length(sCommon)), sCommon, matrix(NA, ncol=6, nrow=length(sCommon))); colnames(sMat) <- colnames(SilRes)
      wMat <- data.frame(rep("Water", length(wCommon)), rep(W_density[j], length(wCommon)),
                         rep(day[k], length(wCommon)), wCommon, matrix(NA, ncol=6, nrow=length(wCommon))); colnames(wMat) <- colnames(SilRes)
      
      SilRes <- rbind(SilRes, rbind(sMat, wMat))
    }
  }
  SilRes_ACR <- SilRes_HH <- SilRes
  
  #-- Silverman test
  for(j in 1:nrow(SilRes)){ #- j=1
    #- Vector of the taxon to be analyzed
    if(SilRes$Source[j] == "Soil"){
      ldf <- sTax[[i]][[which(S_density == SilRes$Density[j])]][[which(day == SilRes$Day[j])]]
    }else{
      ldf <- wTax[[i]][[which(W_density == SilRes$Density[j])]][[which(day == SilRes$Day[j])]]
    }
    ldf <- ldf[,7:ncol(ldf)]
    ldf <- clr(ldf)
    
    ldf <- as.vector(ldf[,colnames(ldf) == SilRes$Taxa[j]])
    SilRes$Count[j] <- sum(ldf>0)/length(ldf); SilRes$Mean[j] <- mean(ldf)
    
    ldf <- robust_scaling(ldf) #- scaling
    
    #- ACR
    Siltest <- modetest(data=ldf, method="ACR", B=5000)
    SilRes_ACR$pvalue[j] <- Siltest$p.value
    SilRes_ACR$Statistic[j] <- Siltest$statistic
    SilRes_ACR$Sta_SD[j] <- Siltest$statistic/sd(ldf)
    SilRes_ACR$Sta_MAD[j] <- Siltest$statistic/mad(ldf)
    
    #- HH
    Siltest <- modetest(data=ldf, method="HH", B=5000)
    SilRes_HH$pvalue[j] <- Siltest$p.value
    SilRes_HH$Statistic[j] <- Siltest$statistic
    SilRes_HH$Sta_SD[j] <- Siltest$statistic/sd(ldf)
    SilRes_HH$Sta_MAD[j] <- Siltest$statistic/mad(ldf)
  }
  
  Pgroup <- data.frame(ACR_p=SilRes_ACR$pvalue, HH_p=SilRes_HH$pvalue)
  ACRq <- (p.adjust(SilRes_ACR$pvalue, method="fdr") < 0.05)
  HHq <- (p.adjust(SilRes_HH$pvalue, method="fdr") < 0.05)
  QQ <- as.factor(ACRq + HHq)
  df <- data.frame(Pgroup, QQ=QQ)

  out1 <- ggplot(df, aes(x=ACR_p, y=HH_p, color=QQ))+
    geom_point(size=0.4)+theme_bw()+
    xlab("p-value (ACR method)")+ylab("p-value (HH method)")+
    geom_abline(slope=1, intercept=0, linetype="dashed", color="blue")+
    theme(text=element_text(size=10))+
    scale_color_manual(values=c("grey70", "orange", "#e31a1c"),
                       labels=c("Both Non-significant", "ACR or HH only", "Both significant"))
  
  df <- data.frame(ACR=SilRes_ACR$Statistic, HH=SilRes_HH$Statistic)
  
  out2 <- ggplot(df, aes(x=ACR, y=HH))+
    geom_point(size=1)+theme_bw()+
    xlab("Multimodality before standard (ACR method)")+ylab("Multimodality before standard (HH method)")+
    theme(text=element_text(size=10))
  
  SilRess[[i]] <- SilRes_HH
}

if(length(list.dirs(sprintf("%s/FigS5", fd), recursive=F))==0){
  dir.create(sprintf("%s/FigS5", fd))
}

pdf(sprintf("%s/FigS5/HHvsACR_CLR_pvalue.pdf", fd), width=4.5, height=3.5)
out1
dev.off()
pdf(sprintf("%s/FigS5/HHvsACR_CLR.pdf", fd), width=3.5, height=3.5)
out2
dev.off()

