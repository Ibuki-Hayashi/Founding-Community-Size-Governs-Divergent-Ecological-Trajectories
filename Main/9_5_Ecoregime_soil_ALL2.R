library(ggplot2); library(vegan);library(gridExtra);library(lemon);library(tidyr);library(markdown)
library(stringr);library(ggforce);library(cowplot);library(ggnewscale); library(reshape2);library(ggtext)
library(ecoregime); library(smacof)
source("functions/functions.R")

# -- read data from previous dir
fd <- "Ecoregime_ALL2"

folder <- c(list.files("Output"), "hoge")
if(any(!(folder %in% fd))){
  dir.create(sprintf("Output/%s", fd))
}
fd <- sprintf("Output/%s", fd)

dd1 <- "Data"; dd2 <- "Data2"

exp <- read.csv(sprintf("%s/expdata_rep.csv", dd1))
data <- readRDS(sprintf("%s/List_data_rep.rds", dd2))

tax <- c("cla","ord","fam","gen","ASV","OTU99","OTU97")
S_density <- c(1, 0.1, 0.01, 0.001)

Inocol <- c("#3B0072", "#903090", "#E735A5", "#F68FF6")
Inocol <- adjustcolor(Inocol, alpha.f = 0.6)

#### -- Data -- ####
#- ndata[[i]][[j]]:
#- i = the number of taxonomy levels, j = the number of Medium

GGG <- function(data, M=S_density){ #- data <- data[[2]]
  data <- data[,!(colnames(data) %in% c("Sample_name", "Plate", "Plate_position"))]
  data <- cbind(data, Random=runif(nrow(data), min=0, max=0.3))
  data <- data[(data$Source=="soil"), ]
  data$Day <- (data$Day)/2
  data$Location <- as.integer(factor(data$Location, levels=unique(data$Location)))
  
  return(data)
}

ndata <- list()
for(i in 1:length(tax)){
  ndata[[i]] <- GGG(data=data[[i]])
}

###--- RETRA-EDR function
#- Minsegs is a hyper parameter
Msegs <- c(5:7) #-- Calculation in multiple Minsegs

for(m in 1:length(Msegs)){ #- m=1
  folder <- c(list.files("Output/Ecoregime_ALL2/figures"), "hoge")
  if(any(!(folder %in% sprintf("Minsegs_%s", Msegs[m])))){
    dir.create(sprintf("Output/Ecoregime_ALL2/figures/Minsegs_%s", Msegs[m]))
  }
  opath <- sprintf("Output/Ecoregime_ALL2/figures/Minsegs_%s", Msegs[m])
  
  #- EDR calculation
  EDR_list <- list(); t_list <- list()
  s_list <- list(); d_list <- list()
  for(i in 6){ #- i=6
    EDR_list[[i]] <- list(); t_list[[i]] <- list()
    s_list[[i]] <- list(); d_list[[i]] <- list()
    
    use <- ndata[[i]]
    ncolor <- c(sum(use$Inoculum_density==S_density[1])/4, sum(use$Inoculum_density==S_density[2])/4,
                sum(use$Inoculum_density==S_density[3])/4, sum(use$Inoculum_density==S_density[4])/4)
    d <- d_list[[i]] <- vegan::vegdist(use[, -c(1:4)])
    states <- s_list[[i]] <- as.integer(use$Day)
    trajectories <- t_list[[i]] <- use$Location
    
    x <- smacofSym(d, ndim=2)$conf
    
    pdf(sprintf("Output/Ecoregime_ALL2/figures/Minsegs_%s/Trajectories_%s_All.pdf", Msegs[m], tax[i]), width=8, height=8)
    plot_edr(x = x, trajectories = trajectories,
             states = states,
             traj.colors = c(rep(Inocol[1], ncolor[1]), rep(Inocol[2], ncolor[2]), rep(Inocol[3], ncolor[3]), rep(Inocol[4], ncolor[4])),
             main = "type = 'trajectories'")
    dev.off()
  }
}
