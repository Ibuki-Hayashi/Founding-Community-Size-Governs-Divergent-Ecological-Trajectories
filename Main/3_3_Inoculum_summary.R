rm(list = ls(all.names = TRUE))

library(ggplot2); library(vegan);library(gridExtra); library(cowplot);library(lemon); library(cluster)
library(tidyr); library(stringr);library(ggforce); library(ggstar);library(ggnewscale); library(ggh4x)

source("functions/functions.R")

#-- information of directory
fd <- "Inoculum"

folder <- c(list.files("Output"), "hoge")
if(any(!(folder %in% fd))){
  dir.create(sprintf("Output/%s", fd))
}
fd <- sprintf("Output/%s", fd)

dd1 <- "Data"; dd2 <- "Data2"

exp <- read.csv(sprintf("%s/expdata_rep.csv", dd1))
S_data <- readRDS(sprintf("%s/Inoculum_data_soil.rds", dd2))
W_data <- readRDS(sprintf("%s/Inoculum_data_water.rds", dd2))

tax <- c("cla","ord","fam","gen","ASV","OTU99","OTU97")
Cli <- c(NA, NA, 10, 15, 20, NA, NA)

Tolong <- function(data, Category, cl){ #- data <- S_data[[1]][[5]]
  if(!(is.data.frame(data))){
    data <- as.data.frame(data)
    if(nrow(data) != 1){
      data <- t(data)
    }
  }
  
  data <- sum_other(data, cl) #- cl=15
  
  ans <- as.data.frame(matrix(NA, ncol=3, nrow=ncol(data)))
  colnames(ans) <- c("Sample_name", "Taxonomy", "Abundance")
  
  ans[,1] <- rep(Category, ncol(data))
  ans[,2] <- colnames(data)
  ans[,3] <- data[1,]
  
  return(ans)
}

#-- 
Long_list <- list(); Long_list[[1]] <- list(); Long_list[[2]] <- list()
Col_pal <- list()

for(i in 3:5){ #- i=5
  Long_list[[1]][[i]] <- s <- Tolong(data=S_data[[1]][[i]], Category="Soil", cl=Cli[i])
  Long_list[[2]][[i]] <- w <-  Tolong(data=W_data[[1]][[i]], Category="Water", cl=Cli[i])
  
  hoge <- rbind(s, w)
  hoge <- hoge[!(hoge[,2] %in% c("Unidentified", "Others")),]
  
  rk <- c()
  for(j in 1:length(unique(hoge[,2]))){ #- i=3
    if(length(hoge[(hoge[,2] == unique(hoge[,2])[j]), 3]) > 1){
      rk[j] <- sum(hoge[(hoge[,2] == unique(hoge[,2])[j]), 3])
    }else{
      rk[j] <- hoge[(hoge[,2] == unique(hoge[,2])[j]), 3]
    }
  }
  
  Col <- unique(hoge[,2]); Ft <- Col[order(rk, decreasing=T)]
  Ft[length(Ft)+1] <- "Unidentified"; Ft[length(Ft)+1] <- "Others"
  
  Col2 <- palettes(Col)
  Col2[length(Col2)+1] <- "#2a333c"; Col2[length(Col2)+1] <- "#d3d3d3"
  
  Col_pal[[i]] <- Col2
  
  Long_list[[1]][[i]][,2] <- hoge <-  factor(Long_list[[1]][[i]][,2], levels=Ft)
  Long_list[[2]][[i]][,2] <- factor(Long_list[[2]][[i]][,2], levels=Ft)
}

#-- plot with Day
Ino_plot <- function(L_i, cal){ #- L_i <- Long_list[[1]][[5]]; cal <- Col_pal[[5]]
  ipal <- cal[levels(L_i$Taxonomy) %in% unique(L_i$Taxonomy)]
  
  ans <- ggplot(L_i)+
    geom_bar(aes(x=Sample_name, y=Abundance, fill=Taxonomy),position="fill",color="black",sta="identity",width=0.95)+
    ylab("Relative abundance")+scale_fill_manual(values=ipal)+
    theme(text=element_text(size=6.4),panel.grid.minor=element_blank(),panel.grid.major=element_blank(),panel.background=element_blank(), axis.text.x=element_blank())+
    scale_y_continuous(expand=c(0,0))+scale_x_discrete(expand=c(0,0))
  
  return(ans)
}

#-- output
S_plot <- list(); W_plot <- list()
for(i in 3:5){
  S_plot[[i]] <- Ino_plot(L_i = Long_list[[1]][[i]], cal <- Col_pal[[i]])
  W_plot[[i]] <- Ino_plot(L_i = Long_list[[2]][[i]], cal <- Col_pal[[i]])
}

#- Soil <- plot_grid(plot_grid(S_plot[[3]], S_plot[[5]], nrow=1), S_plot[[4]], nrow=2)
Soil <- plot_grid(S_plot[[3]], S_plot[[4]] ,S_plot[[5]], nrow=1, rel_widths = c(1.1, 1.5, 1.15), labels=c("Family", "Genus", "ASV"), label_x=1)
Water <- plot_grid(W_plot[[3]], W_plot[[4]] ,W_plot[[5]], nrow=1, rel_widths = c(1.15, 1.82, 1.25), labels=c("Family", "Genus", "ASV"), label_x=1)

Ans <- plot_grid(Soil, NA, Water, nrow=3, rel_heights = c(1,0.05,1))

pdf(sprintf("%s/Inoculum_bar.pdf", fd), width=7, height=8.5)
Ans
dev.off()