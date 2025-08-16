rm(list = ls(all.names = TRUE))

library(ggplot2); library(vegan);library(gridExtra); library(cowplot);library(lemon); library(cluster)
library(tidyr); library(stringr);library(ggforce); library(ggstar);library(ggnewscale); library(ggh4x)
library(multimode); library(MCMCpack); library(doParallel); library(ggplate)

source("Functions/functions.R")

#-- information of directory
fd <- "Well_position"

folder <- c(list.files("Output"), "hoge")
if(any(!(folder %in% fd))){
  dir.create(sprintf("Output/%s", fd))
}
fd <- sprintf("Output/%s", fd)

dd1 <- "Data"; dd2 <- "Data2"

exp <- read.csv(sprintf("%s/expdata_rep.csv", dd1))
data <- readRDS(sprintf("%s/List_data_rep.rds", dd2))
tax <- c("cla","ord","fam","gen","ASV","OTU99","OTU97")
Density <- c(1, 0.1, 0.01, 0.001)
day <- c(2, 4, 6, 8)
scal <- 1.5

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

#-- Clustering
Figin <- function(df, scaling=scal){ #- df = sTax[[6]][[2]][[4]]
  mat <- df[,7:ncol(df)]
  ks <- ksil(mat)
  
  dist <- vegdist(mat, method="bray")
  pcoa <- cmdscale(dist, eig = TRUE, k = 2)  # PCoA
  
  outmat <- data.frame(Sample = rownames(mat),
                             Axis1 = pcoa$points[,1],
                             Axis2 = pcoa$points[,2],
                             Cluster = factor(ks$clustering, levels=sort(unique(ks$clustering))))
  
  pcplot <- ggplot(outmat, aes(x=Axis1, y=Axis2, color=Cluster))+
    geom_point(size=scaling)+
    theme_classic()+
    scale_color_manual(values=palette_30[1:length(unique(ks$clustering))])+
    theme(text = element_text(size=7*scaling))+
    ggtitle(sprintf("%s-Density:%s Day:%s", df$Source[1], 
                    df$Inoculum_density[1], 
                    df$Day[1]))
  
  ndf <- data.frame(Location=df$Location,
                    Cluster=as.character(factor(ks$clustering, levels=sort(unique(ks$clustering)))))
  tndf <- as_tibble(ndf)
  tndf$Location <- sapply(as.character(tndf$Location),
                          function(value){
                            v1 <- toupper(str_sub(value,1,1))
                            return(paste(v1,str_sub(value,2,nchar(value)),sep=""))})
  
  plate <- plate_plot(
    data=tndf,
    position=Location,
    plate_size=384,
    value=Cluster,
    colour=palette_30[1:length(unique(ks$clustering))],
    title_size=8*scaling)+
    theme(legend.text = element_text(size=7*scaling),
          legend.title = element_text(size=6*scaling),
          axis.text = element_text(size=6*scaling),
          plot.title = element_text(size=8*scaling))
  
  blank <- ggplot() + theme_void()
  
  ans <- plot_grid(pcplot, plate, rel_widths = c(0.7, 1), ncol=2)
}

Spl <- list(); Wpl <- list()
for(j in 1:length(Density)){
  Spl[[j]] <- list(); Wpl[[j]] <- list()
  for(k in 1:length(day)){
    Spl[[j]][[k]] <- Figin(df=sTax[[6]][[j]][[k]])
    Wpl[[j]][[k]] <- Figin(df=wTax[[6]][[j]][[k]])
  }
}
blank <- ggplot() + theme_void()

S_day8 <- plot_grid(Spl[[1]][[4]], blank, Spl[[2]][[4]], blank,
                    Spl[[3]][[4]], blank, Spl[[4]][[4]],
                    ncol=1, nrow=8, rel_heights = c(1, 0.1, 1, 0.1, 1, 0.1, 1, 0.1))
W_day8 <- plot_grid(Wpl[[1]][[4]], blank, Wpl[[2]][[4]], blank,
                    Wpl[[3]][[4]], blank, Wpl[[4]][[4]],
                    ncol=1, nrow=8, rel_heights = c(1, 0.1, 1, 0.1, 1, 0.1, 1, 0.1))

ggsave(S_day8, filename=sprintf("%s/OTU99_Soil_cluster.pdf", fd), height=1100*scal/5, width=800*scal/5, units="mm")
ggsave(W_day8, filename=sprintf("%s/OTU99_Water_cluster.pdf", fd), height=1100*scal/5, width=800*scal/5, units="mm")
