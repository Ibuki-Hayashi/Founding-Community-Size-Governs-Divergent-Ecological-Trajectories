library(ggplot2); library(vegan);library(gridExtra); library(cowplot);library(lemon)
library(tidyr); library(stringr);library(ggforce); library(ggnewscale)
library(ggtext)

#-- information of directory
fd <- "Alpha"

folder <- c(list.files("Output"), "hoge")
if(any(!(folder %in% fd))){
  dir.create(sprintf("Output/%s", fd))
}
fd <- sprintf("Output/%s", fd)

exp <- read.csv("Data/expdata_rep.csv")

#-- read data from previous dir
data <- readRDS("Output/List_data_rep.rds"); tax <- c("Class","Order","Family","Genus","ASV","OTU99","OTU97")
ti <- c(2,4,6,8)
Daycol <- c("#1f78b4", "#33a02c", "#ff7f00", "#6a3d9a")

# -- Richness
Richness_01 <- function(data, thr=0){
  ans1 <- c()
  for(i in 1:nrow(data)){
    ans1[i] <- 0
    for(j in 1:ncol(data)){
      if(data[i,j] != thr){(ans1[i] <- ans1[i]+1)}
    }
  }
  return(ans1)
}

#-- calcurate alphas
alpha_mat <- function(data, enumber, ino){
  seqtab <- data[(data$Source == ino),(enumber+1):ncol(data)]
  exp <- data[(data$Source == ino), 1:enumber]
  
  ans <- data.frame(Shannon = exp(diversity(seqtab, index="shannon")),
                    Richness = Richness_01(seqtab))
  ans <- cbind(exp, ans); return(ans)
}

almat_ls_s <- list(); almat_ls_w <- list()
for(i in 3:7){
  almat_ls_s[[i]] <- alpha_mat(data=data[[i]], enumber=6, ino="soil")
  almat_ls_s[[i]]$Inoculum_density <- factor(almat_ls_s[[i]]$Inoculum_density, levels=c(1, 0.1, 0.01, 0.001))
  almat_ls_s[[i]]$Day <- as.factor(almat_ls_s[[i]]$Day)
  almat_ls_w[[i]] <- alpha_mat(data=data[[i]], enumber=6, ino="water")
  almat_ls_w[[i]]$Inoculum_density <- factor(almat_ls_w[[i]]$Inoculum_density, levels=c(1, 0.1, 0.01, 0.001))
  almat_ls_w[[i]]$Day <- as.factor(almat_ls_w[[i]]$Day)
}

# -- violinplot
#box_med
box_rich <- function(ddd, title){
  box <- ggplot(ddd, aes(x=Day, y=Richness))+
    geom_violin(aes(group=Day))+ylab("Number of taxa of community<br>(Species richness)")+
    geom_jitter(aes(color=Day), size=1)+
    theme_bw()+theme(text=element_text(size=10),axis.title.y=element_markdown(),
                     legend.position="none",
                     strip.background = element_rect(fill = "black"),
                     strip.text = element_text(color = "white", face = "bold"))+
    scale_color_manual(values=Daycol)+
    facet_wrap(~Inoculum_density, nrow=2, ncol=2)+ggtitle(title)
  return(box)
}

box_shannon <- function(ddd, title){
  box <- ggplot(ddd, aes(x=Day, y=Shannon))+
    geom_violin(aes(group=Day))+ylab("<i>&alpha;</i> diversity of community<br>(Shannon effective number)")+
    geom_jitter(aes(color=Day), size=1)+
    theme_bw()+theme(text=element_text(size=10),axis.title.y=element_markdown(),
                     legend.position="none",
                     strip.background = element_rect(fill = "black"),
                     strip.text = element_text(color = "white", face = "bold"))+
    scale_color_manual(values=Daycol)+
    facet_wrap(~Inoculum_density, nrow=2, ncol=2)+ggtitle(title)
  return(box)
}

ssh <- list(); sri <- list()
wsh <- list(); wri <- list()
su <- list(); wu <- list(); au <- list()
for(i in 3:7){
  ssh[[i]] <- box_shannon(almat_ls_s[[i]], sprintf("%s_Soil", tax[i]))
  sri[[i]] <- box_rich(almat_ls_s[[i]], sprintf("%s_Soil", tax[i]))
  su[[i]] <- plot_grid(ssh[[i]], sri[[i]], ncol=2)
  
  wsh[[i]] <- box_shannon(almat_ls_w[[i]], sprintf("%s_Water", tax[i]))
  wri[[i]] <- box_rich(almat_ls_w[[i]], sprintf("%s_Water", tax[i]))
  wu[[i]] <- plot_grid(wsh[[i]], wri[[i]], ncol=2)
  
  au[[i]] <- plot_grid(su[[i]], wu[[i]], nrow=2)
}

# -- output
pdf(sprintf("%s/Soil_alphas_384_effective.pdf", fd),width=8,height=6)
su[[3]]; su[[4]]; su[[5]]; su[[6]]; su[[7]]
dev.off()

pdf(sprintf("%s/Water_alphas_384_effective.pdf", fd),width=8,height=6)
wu[[3]]; wu[[4]]; wu[[5]]; wu[[6]]; wu[[7]]
dev.off()

pdf(sprintf("%s/All_alphas_384_effective.pdf", fd),width=8,height=7)
au[[3]]; au[[4]]; au[[5]]; au[[6]]; au[[7]]
dev.off()

pdf(sprintf("%s/ForFig_alphas_384_effective.pdf", fd),width=8,height=7)
au[[6]]
dev.off()