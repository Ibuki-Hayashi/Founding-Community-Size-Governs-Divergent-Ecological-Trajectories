source("Functions/functions.R")

#-- information of directory
fd <- "Output"; dd <- "Data"

#-- information of directory
exp <- read.csv(sprintf("%s/expdata_rep.csv", dd))
data <- readRDS(sprintf("%s/List_data_rep.rds", fd))
Sil <- readRDS(sprintf("%s/Binomialfir_Master_ACR.rds", fd))

tax <- c("cla","ord","fam","gen","ASV","OTU99","OTU97")
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

#-- Silplot
corpl1 <- list(); corpl2 <- list()
for(i in 6){ #- i=6
  SilRes <- Sil[[i]]
  
  uSil <- SilRes[(SilRes$LL != Inf),]
  uSil$Act_CV <- uSil$sigma/uSil$Mean
  uSil$Sta_M <- uSil$Statistic/uSil$Mean
  uSil$Ino_CV <- log(uSil$Ino_CV, 10)
  uSil$unique <- sprintf("%s_%s", uSil$Source, uSil$Taxa)
  
  corpl1[[i]] <- gList(); corpl2[[i]] <- gList()
  itr <- 1
  for(j in 1:length(unique(uSil$unique))){ #- j=1
    df <- uSil[uSil$unique == unique(uSil$unique)[j],]
    key <- c()
    for(k in 1:4){
      key[k] <- length(df$Day[df$Day == day[k]])
    }
    if(sum(key > 2) >= 2){
      corpl1[[i]][[itr]] <- ggplot(df, aes(x=Ino_CV, y=CV, color=as.factor(Day)))+
        geom_point(size=1.3)+geom_line(size=0.25)+theme_bw()+xlab("")+ylab("")+
        ggtitle(sprintf("%s", unique(uSil$unique)[j]))+theme(legend.position="none", text=element_text(size=7.2))
      corpl2[[i]][[itr]] <- ggplot(df, aes(x=Ino_CV, y=Statistic, color=as.factor(Day)))+
        geom_point(size=1.3)+geom_line(size=0.25)+theme_bw()+xlab("")+ylab("")+
        ggtitle(sprintf("%s", unique(uSil$unique)[j]))+theme(legend.position="none", text=element_text(size=7.2))
      
      itr <- itr+1
    }
  }
  leg <- ggplot(data=data.frame(x=c(2,4,6,8), Day=as.factor(c(2,4,6,8))))+geom_point(aes(x=x, y=Day, color=Day))+theme_bw()
  leg <- g_legend(leg)
}

n991 <- grid.arrange(grobs=corpl1[[6]], ncol=6, nrow=2)
n991 <- plot_grid(n991, leg, ncol=2, rel_widths=c(1, 0.1))
n992 <- grid.arrange(grobs=corpl2[[6]], ncol=6, nrow=2)
n992 <- plot_grid(n992, leg, ncol=2, rel_widths=c(1, 0.1))
n99 <- plot_grid(NA, n991, n992, ncol=1, rel_heights=c(0.05, 1, 1), labels=c("OTU99", NA, NA))

#-- Save
if(length(list.dirs(sprintf("%s/FigS16", fd), recursive=F))==0){
  dir.create(sprintf("%s/FigS16", fd))
}

pdf(sprintf("%s/FigS16/EachTaxa_Cor.pdf", fd), width=8.2, height=4)
plot(n99)
dev.off()
