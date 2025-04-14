library(ggplot2); library(vegan);library(gridExtra);library(lemon);library(tidyr);library(markdown)
library(stringr);library(ggforce);library(cowplot);library(ggnewscale); library(reshape2);library(ggtext)
library(ecoregime)
source("Functions/functions.R")

# -- read data from previous dir
dd <- "Data"; fd <- "Output"
dir.create(sprintf("%s/Ecoregime_soil", fd))

exp <- read.csv(sprintf("%s/expdata_rep.csv"), dd)
data <- readRDS(sprintf("%s/List_data_rep.rds", fd))

tax <- c("cla","ord","fam","gen","ASV","OTU99","OTU97")
S_density <- c(1, 0.1, 0.01, 0.001)

#### -- Data -- ####
#- ndata[[i]][[j]]:
#- i = the number of taxonomy levels, j = the number of Medium

GGG <- function(data, M=S_density){ #- data <- data[[2]]
  data <- data[,!(colnames(data) %in% c("Sample_name", "Plate", "Plate_position"))]
  data <- cbind(data, Random=runif(nrow(data), min=0, max=0.3))
  data <- data[(data$Source=="soil"), ]
  data$Day <- (data$Day)/2
  
  ans <- list()
  for(i in 1:length(M)){ #- i=1
    ans[[i]] <- data[(data$Inoculum_density == M[i]),]
    ans[[i]]$Location <- as.integer(factor(ans[[i]]$Location, levels=unique(ans[[i]]$Location)))
  }
  
  return(ans)
}

ndata <- list()
for(i in 1:length(tax)){
  ndata[[i]] <- GGG(data=data[[i]])
}

###--- RETRA-EDR function
#- Minsegs is a hyper parameter
Msegs <- c(5:7) #-- Calculation in multiple Minsegs

for(m in 1:length(Msegs)){ #- m=1
  folder <- c(list.files(sprintf("%s/Ecoregime_soil", fd)), "hoge")
  if(any(!(folder %in% sprintf("Minsegs_%s", Msegs[m])))){
    dir.create(sprintf("%s/Ecoregime_soil/Minsegs_%s", fd, Msegs[m]))
  }
  opath <- sprintf("%s/Ecoregime_soil/Minsegs_%s", fd, Msegs[m])
  
  #- EDR calculation
  EDR_list <- list(); t_list <- list()
  s_list <- list(); d_list <- list()
  for(i in 4:length(tax)){ #- i=5
    EDR_list[[i]] <- list(); t_list[[i]] <- list()
    s_list[[i]] <- list(); d_list[[i]] <- list()
    for(j in 1:length(S_density)){ #- j=1
      #- calculation part
      use <- ndata[[i]][[j]]
      
      d <- d_list[[i]][[j]] <- vegan::vegdist(use[, -c(1:4)])
      trajectories <- t_list[[i]][[j]] <- use$Location
      states <- s_list[[i]][[j]] <- as.integer(use$Day)
      
      EDR_list[[i]][[j]] <- retra_edr(d=d, trajectories=trajectories, states=states, minSegs=Msegs[m])
    }
  }
  
  #- EDR trajectories plots
  for(i in 4:length(tax)){
    pdf(sprintf("%s/Ecoregime_soil/Minsegs_%s/Trajectories_%s.pdf", fd, Msegs[m], tax[i]), width=16, height=4.4)
    par(mfrow = c(1, 4))
    for(j in 1:length(S_density)){
      srn <- 1; srt <- "T1"
      for(t in 2:length(EDR_list[[i]][[j]])){
        if(EDR_list[[i]][[j]][[srn]]$Length < EDR_list[[i]][[j]][[t]]$Length){
          srn <- t; srt <- sprintf("T%s",as.character(t))
        }
      }
      plot(x = EDR_list[[i]][[j]], d = d_list[[i]][[j]], trajectories = t_list[[i]][[j]], states = s_list[[i]][[j]],
           select_RT = srt, #if we need to highlight some RT
           traj.colors = "lightblue", RT.colors = "orange", sel.color = "darkgreen",
           link.lty = 1, asp = 1, main = sprintf("Dilution rate : %s", S_density[j]))
    }
    dev.off()
    
    #- EDR dynamic metrices
    ans <- matrix(NA, ncol=length(S_density), nrow=3)
    for(j in 1:length(S_density)){
      rownames(ans) <- c("dDis","dBD","dEv"); colnames(ans) <- paste("Density", S_density, sep="_")
      srn <- 1; srt <- "T1"
      for(t in 2:length(EDR_list[[i]][[j]])){
        if(EDR_list[[i]][[j]][[srn]]$Length < EDR_list[[i]][[j]][[t]]$Length){
          srn <- t; srt <- sprintf("T%s",as.character(t))
        }
      }
      d = d_list[[i]][[j]]; trajectories = t_list[[i]][[j]]; states = s_list[[i]][[j]]
      ans[1,j] <- dDis(d = d, d.type = "dStates", trajectories = trajectories, states = states, reference = as.character(srn))
      ans[2,j] <- dBD(d = d, d.type = "dStates", trajectories = trajectories, states = states)
      ans[3,j] <- dEve(d = d, d.type = "dStates", trajectories = trajectories, states = states)
    }
    write.csv(ans, sprintf("%s/Ecoregime_soil/Minsegs_%s/Metrices_%s.csv", fd, Msegs[m], tax[i]))
  }
}
