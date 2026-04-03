source("Functions/functions.R")
#-- information of directory
fd <- "Output"; dd <- "Data"

exp <- read.csv(sprintf("%s/expdata_rep.csv", dd))
data <- readRDS(sprintf("%s/List_data_rep.rds", fd)); tax <- c("cla","ord","fam","gen","ASV","OTU99","OTU97")
S_density <- c(1, 0.1, 0.01, 0.001)
W_density <- c(1, 0.1, 0.01, 0.001)

#-- make data ordered, & long data & color palettes
thr1 <- c(8,12,14,18,36,22,28)
thr2 <- c(10,16,20,24,48,30,42)

ans <- list()
width_fig <- c(0.15, 0.15, 0.2, 0.4, 0.2, 0.15, 0.15)
legncol <- c(1,1,1,1,2,1,1)
for(i in 6){ #- i=6
  df <- data[[i]]
  df$Day <- sprintf("Day%s", df$Day)
  dfs <- df[df$Source=="soil",]; dfw <- df[df$Source=="water",]
  s_key <- colSums(dfs[,7:ncol(dfs)]); w_key <- colSums(dfw[,7:ncol(dfw)]); a_key <- colSums(df[,7:ncol(df)])
  
  s_tax <- names(s_key)[s_key > sort(s_key, TRUE)[thr1[i]/2]-1]
  w_tax <- names(w_key)[w_key > sort(w_key, TRUE)[thr1[i]/2]-1]
  
  a_key <- a_key[!(names(a_key) %in% union(s_tax, w_tax))]
  a_tax <- names(a_key)[a_key > sort(a_key, TRUE)[thr2[i]-thr1[i]]-1]
  
  u_tax <- unique(c(s_tax, w_tax, a_tax))
  u_tax <- u_tax[!(u_tax %in% c("Unidentified", "Others"))]
  u_tax <- c(u_tax, "Unidentified", "Others")
  cpal <- palettes(u_tax[-c(1,2)]); cpal <- c(cpal, "#2a333c", "#d3d3d3")
  
  #- Make others
  mt <- df[,7:ncol(df)]
  ndf <- cbind(df[,1:6], cbind(mt[, colnames(mt) %in% u_tax], Others=rowSums(mt[, !(colnames(mt) %in% u_tax)])))
  lndf <- gather(ndf, Taxonomy, Abundance, -c(1:6))
  lndf$Taxonomy <- factor(lndf$Taxonomy, levels=u_tax)
  
  #- Location is the order
  dfs_well <- dfs[(dfs$Day=="Day8"), c(7:ncol(dfs))] %>%
    vegdist(method="bray") %>% hclust(method="average")
  sk <- dfs$Location[dfs_well[["order"]]]
  
  dfw_well <- dfw[(dfw$Day=="Day8"), c(7:ncol(dfw))] %>%
    vegdist(method="bray") %>% hclust(method="average")
  wk <- dfw$Location[dfw_well[["order"]]]
  
  s_plot <- list(); w_plot <- list()
  for(j in 1:length(S_density)){ #- j=1
    wdf <- lndf[(lndf$Source=="water")&(lndf$Inoculum_density==W_density[j]),]
    wdf$Location <- factor(wdf$Location, levels=wk)
    wdf <- wdf[wdf$Day == "Day2",]
    
    w_plot[[j]] <- ggplot(wdf)+
      geom_bar(aes(x=Location, y=Abundance, fill=Taxonomy), position="fill", sta="identity", width=0.9)+
      ggtitle(sprintf("Water_%s_Day2", S_density[j]))+xlab("Replicate communities")+ylab("Relative abundance")+
      theme(text = element_text(size=7), panel.grid.minor=element_blank(), axis.ticks=element_blank(),
            panel.spacing=unit(1.5,"pt"), plot.margin=unit(c(1,2.5,1,2.5),"pt"),
            panel.grid.major=element_blank(), panel.background=element_blank(), axis.text.x=element_blank(),strip.text.x=element_text(margin=margin(0,0,0,0,unit="pt")))+
      scale_y_continuous(expand=c(0,0))+scale_x_discrete(expand=c(0,0))+
      scale_fill_manual(values=cpal)+guides(fill = guide_legend(ncol = legncol[i]))
    w_plot[[j]] <- w_plot[[j]]+theme(legend.position = "none")
  }
  
  ans[[i]] <- grid.arrange(w_plot[[1]], w_plot[[2]], w_plot[[3]], w_plot[[4]], nrow=1, ncol=4)
}

#- output
if(length(list.dirs(sprintf("%s/Fig3", fd), recursive=F))==0){
  dir.create(sprintf("%s/Fig3", fd))
}

pdf(sprintf("%s/Fig3/BarFig3.pdf", fd), height=1, width=8)
plot(ans[[6]])
dev.off()

