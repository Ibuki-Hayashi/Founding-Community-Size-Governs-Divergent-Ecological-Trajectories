library(ggplot2); library(vegan);library(gridExtra);library(lemon);library(tidyr);library(markdown)
library(stringr);library(ggforce);library(cowplot);library(ggnewscale); library(reshape2);library(ggtext)
library(ecoregime)
source("functions/functions.R")

# -- read data from previous dir
fd <- "Ecoregime_water"

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

#### -- Data -- ####
#- ndata[[i]][[j]]:
#- i = the number of taxonomy levels, j = the number of Medium

GGG <- function(data, M=S_density){ #- data <- data[[2]]
  data <- data[,!(colnames(data) %in% c("Sample_name", "Plate", "Plate_position"))]
  data <- cbind(data, Random=runif(nrow(data), min=0, max=0.3))
  data <- data[(data$Source=="water"), ]
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
  folder <- c(list.files("Output/Ecoregime_water/figures"), "hoge")
  if(any(!(folder %in% sprintf("Minsegs_%s", Msegs[m])))){
    dir.create(sprintf("Output/Ecoregime_water/figures/Minsegs_%s", Msegs[m]))
  }
  opath <- sprintf("Output/Ecoregime_water/figures/Minsegs_%s", Msegs[m])
  
  #- EDR calculation
  EDR_list <- list(); t_list <- list()
  s_list <- list(); d_list <- list()
  for(i in 6){ #- i=5
    EDR_list[[i]] <- list(); t_list[[i]] <- list()
    s_list[[i]] <- list(); d_list[[i]] <- list()
    for(j in 1:length(S_density)){ #- j=1
      #- calculation part
      use <- ndata[[i]][[j]]
      
      d <- d_list[[i]][[j]] <- vegan::vegdist(use[, -c(1:4)])
      trajectories <- t_list[[i]][[j]] <- use$Location
      states <- s_list[[i]][[j]] <- as.integer(use$Day)
      
      set.seed(111111)
      EDR_list[[i]][[j]] <- retra_edr(d=d, trajectories=trajectories, states=states, minSegs=Msegs[m])
    }
  }
  
  #- EDR trajectories plots
  for(i in 6){ #- i=6
    ###########- Max size
    #- EDR trajectories
    pdf(sprintf("Output/Ecoregime_water/figures/Minsegs_%s/Trajectories_%s_Size.pdf", Msegs[m], tax[i]), width=16, height=4.4)
    par(mfrow = c(1, 4))
    for(j in 1:length(S_density)){ #- j=2
      summ <- summary(EDR_list[[i]][[j]])
      summ <- summ[(is.na(summ$Size) == FALSE), ] #-- remove NA
      srt <- summ$ID[which.max(summ$Size)] #-- select max length
      srn <- as.numeric(str_sub(srt, 2,-1))
      
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
      
      summ <- summary(EDR_list[[i]][[j]])
      summ <- summ[(is.na(summ$Size) == FALSE), ] #-- remove NA
      srt <- summ$ID[which.max(summ$Size)] #-- select max length
      srn <- as.numeric(str_sub(srt, 2,-1))
      d = d_list[[i]][[j]]; trajectories = t_list[[i]][[j]]; states = s_list[[i]][[j]]
      
      ans[1,j] <- dDis(d = d, d.type = "dStates", trajectories = trajectories, states = states, reference = as.character(srn))
      ans[2,j] <- dBD(d = d, d.type = "dStates", trajectories = trajectories, states = states)
      ans[3,j] <- dEve(d = d, d.type = "dStates", trajectories = trajectories, states = states)
    }
    write.csv(ans, sprintf("Output/Ecoregime_water/figures/Minsegs_%s/Metrices_%s_Size.csv", Msegs[m], tax[i]))
    
    ###########- Min Average Link
    #- EDR trajectories
    pdf(sprintf("Output/Ecoregime_water/figures/Minsegs_%s/Trajectories_%s_Link.pdf", Msegs[m], tax[i]), width=16, height=4.4)
    par(mfrow = c(1, 4))
    for(j in 1:length(S_density)){ #- j=3
      summ <- summary(EDR_list[[i]][[j]])
      summ <- summ[(is.na(summ$Avg_link) == FALSE), ] #-- remove NA
      srt <- summ$ID[which.min(summ$Avg_link)] #-- select max length
      if(length(srt) > 1){
        srt <- srt[1] #-- if there are multiple max, select the first one
      }
      srn <- as.numeric(str_sub(srt, 2,-1))
      
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
      
      summ <- summary(EDR_list[[i]][[j]])
      summ <- summ[(is.na(summ$Avg_link) == FALSE), ] #-- remove NA
      srt <- summ$ID[which.min(summ$Avg_link)] #-- select max length
      if(length(srt) > 1){
        srt <- srt[1] #-- if there are multiple max, select the first one
      }
      srn <- as.numeric(str_sub(srt, 2,-1))
      d = d_list[[i]][[j]]; trajectories = t_list[[i]][[j]]; states = s_list[[i]][[j]]
      
      ans[1,j] <- dDis(d = d, d.type = "dStates", trajectories = trajectories, states = states, reference = as.character(srn))
      ans[2,j] <- dBD(d = d, d.type = "dStates", trajectories = trajectories, states = states)
      ans[3,j] <- dEve(d = d, d.type = "dStates", trajectories = trajectories, states = states)
    }
    write.csv(ans, sprintf("Output/Ecoregime_water/figures/Minsegs_%s/Metrices_%s_Link.csv", Msegs[m], tax[i]))
    
    ###########- Max density
    #- EDR trajectories
    pdf(sprintf("Output/Ecoregime_water/figures/Minsegs_%s/Trajectories_%s_Density.pdf", Msegs[m], tax[i]), width=16, height=4.4)
    par(mfrow = c(1, 4))
    for(j in 1:length(S_density)){
      summ <- summary(EDR_list[[i]][[j]])
      summ <- summ[(is.na(summ$Avg_density) == FALSE), ] #-- remove NA
      srt <- summ$ID[which.max(summ$Avg_density)] #-- select max length
      if(length(srt) > 1){
        srt <- srt[1] #-- if there are multiple max, select the first one
      }
      srn <- as.numeric(str_sub(srt, 2,-1))
      
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
      
      summ <- summary(EDR_list[[i]][[j]])
      summ <- summ[(is.na(summ$Avg_density) == FALSE), ] #-- remove NA
      srt <- summ$ID[which.max(summ$Avg_density)] #-- select max length
      if(length(srt) > 1){
        srt <- srt[1] #-- if there are multiple max, select the first one
      }
      srn <- as.numeric(str_sub(srt, 2,-1))
      d = d_list[[i]][[j]]; trajectories = t_list[[i]][[j]]; states = s_list[[i]][[j]]
      
      ans[1,j] <- dDis(d = d, d.type = "dStates", trajectories = trajectories, states = states, reference = as.character(srn))
      ans[2,j] <- dBD(d = d, d.type = "dStates", trajectories = trajectories, states = states)
      ans[3,j] <- dEve(d = d, d.type = "dStates", trajectories = trajectories, states = states)
    }
    write.csv(ans, sprintf("Output/Ecoregime_water/figures/Minsegs_%s/Metrices_%s_Density.csv", Msegs[m], tax[i]))
    
    ###########- Max depth
    #- EDR trajectories
    pdf(sprintf("Output/Ecoregime_water/figures/Minsegs_%s/Trajectories_%s_Depth.pdf", Msegs[m], tax[i]), width=16, height=4.4)
    par(mfrow = c(1, 4))
    for(j in 1:length(S_density)){
      summ <- summary(EDR_list[[i]][[j]])
      summ <- summ[(is.na(summ$Avg_depth) == FALSE), ] #-- remove NA
      srt <- summ$ID[which.max(summ$Avg_depth)] #-- select max length
      if(length(srt) > 1){
        srt <- srt[1] #-- if there are multiple max, select the first one
      }
      srn <- as.numeric(str_sub(srt, 2,-1))
      
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
      
      summ <- summary(EDR_list[[i]][[j]])
      summ <- summ[(is.na(summ$Avg_depth) == FALSE), ] #-- remove NA
      srt <- summ$ID[which.max(summ$Avg_depth)] #-- select max length
      if(length(srt) > 1){
        srt <- srt[1] #-- if there are multiple max, select the first one
      }
      srn <- as.numeric(str_sub(srt, 2,-1))
      d = d_list[[i]][[j]]; trajectories = t_list[[i]][[j]]; states = s_list[[i]][[j]]
      
      ans[1,j] <- dDis(d = d, d.type = "dStates", trajectories = trajectories, states = states, reference = as.character(srn))
      ans[2,j] <- dBD(d = d, d.type = "dStates", trajectories = trajectories, states = states)
      ans[3,j] <- dEve(d = d, d.type = "dStates", trajectories = trajectories, states = states)
    }
    write.csv(ans, sprintf("Output/Ecoregime_water/figures/Minsegs_%s/Metrices_%s_Depth.csv", Msegs[m], tax[i]))
    
    ###########- Max rank
    #- EDR trajectories
    pdf(sprintf("Output/Ecoregime_water/figures/Minsegs_%s/Trajectories_%s_Rank.pdf", Msegs[m], tax[i]), width=16, height=4.4)
    par(mfrow = c(1, 4))
    for(j in 1:length(S_density)){ #- j=1
      summ <- summary(EDR_list[[i]][[j]])
      summ <- summ[(is.na(summ$Size) == FALSE)&(is.na(summ$Avg_link) == FALSE)&(is.na(summ$Avg_density) == FALSE)&(is.na(summ$Avg_depth) == FALSE), ] #-- remove NA
      
      summ$Rank <- rank(summ$Size, ties.method = "min") + rank(-summ$Avg_link, ties.method = "min")+
        rank(summ$Avg_density, ties.method = "min") + rank(summ$Avg_depth, ties.method = "min")
      
      srt <- summ$ID[which.min(summ$Rank)] #-- select max length
      if(length(srt) > 1){
        srt <- srt[1] #-- if there are multiple max, select the first one
      }
      srn <- as.numeric(str_sub(srt, 2,-1))
      
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
      
      summ <- summary(EDR_list[[i]][[j]])
      summ <- summ[(is.na(summ$Size) == FALSE)&(is.na(summ$Avg_link) == FALSE)&(is.na(summ$Avg_density) == FALSE)&(is.na(summ$Avg_depth) == FALSE), ] #-- remove NA
      
      summ$Rank <- rank(summ$Size, ties.method = "min") + rank(-summ$Avg_link, ties.method = "min")+
        rank(summ$Avg_density, ties.method = "min") + rank(summ$Avg_depth, ties.method = "min")
      
      srt <- summ$ID[which.min(summ$Rank)] #-- select max length
      if(length(srt) > 1){
        srt <- srt[1] #-- if there are multiple max, select the first one
      }
      srn <- as.numeric(str_sub(srt, 2,-1))
      
      d = d_list[[i]][[j]]; trajectories = t_list[[i]][[j]]; states = s_list[[i]][[j]]
      
      ans[1,j] <- dDis(d = d, d.type = "dStates", trajectories = trajectories, states = states, reference = as.character(srn))
      ans[2,j] <- dBD(d = d, d.type = "dStates", trajectories = trajectories, states = states)
      ans[3,j] <- dEve(d = d, d.type = "dStates", trajectories = trajectories, states = states)
    }
    write.csv(ans, sprintf("Output/Ecoregime_water/figures/Minsegs_%s/Metrices_%s_Rank.csv", Msegs[m], tax[i]))
  }
}
