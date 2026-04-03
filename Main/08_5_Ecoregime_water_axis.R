source("Functions/functions.R")

#-- information of directory
fd <- "Output"; dd <- "Data"

exp <- read.csv(sprintf("%s/expdata_rep.csv", dd))
data <- readRDS(sprintf("%s/List_data_rep.rds", fd))
tax <- c("cla","ord","fam","gen","ASV","OTU99","OTU97")
W_density <- S_density <- c(1, 0.1, 0.01, 0.001)


#### -- Data -- ####
#- ndata[[i]][[j]]:
#- i = the number of taxonomy levels, j = the number of Medium

GGG <- function(data, M=W_density){ #- data <- data[[2]]
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
Msegs <- c(5) #-- Calculation in multiple Minsegs

if(length(list.dirs(sprintf("%s/FigS20water", fd), recursive=F))==0){
  dir.create(sprintf("%s/FigS20water", fd))
}

for(m in 1:length(Msegs)){ #- m=1
  folder <- c(list.files("Output/FigS20water"), "hoge")
  if(any(!(folder %in% sprintf("Minsegs_%s", Msegs[m])))){
    dir.create(sprintf("Output/FigS20water/Minsegs_%s", Msegs[m]))
  }
  opath <- sprintf("Output/FigS20water/Minsegs_%s", Msegs[m])
  
  #- EDR calculation
  EDR_list <- list(); t_list <- list()
  s_list <- list(); d_list <- list()
  for(i in 6){ #- i=4
    EDR_list[[i]] <- list(); t_list[[i]] <- list()
    s_list[[i]] <- list(); d_list[[i]] <- list()
    for(j in 1:length(W_density)){ #- j=2
      #- calculation part
      use <- ndata[[i]][[j]]
      
      d <- d_list[[i]][[j]] <- vegan::vegdist(use[, -c(1:4)])
      
      trajectories <- t_list[[i]][[j]] <- use$Location
      states <- s_list[[i]][[j]] <- as.integer(use$Day)
      
      EDR_list[[i]][[j]] <- retra_edr(d=d, trajectories=trajectories, states=states, minSegs=Msegs[m])
    }
  }
  
  #- EDR trajectories plots
  for(i in 6){ #- i=6
    hdlist <- lapply(d_list[[i]], smacofSym)
    poi1 <- lapply(hdlist, function(x) c(x$conf[,1])); poi1a <- unlist(poi1)
    xMa <- max(poi1a); xMi <- min(poi1a)
    poi2 <- lapply(hdlist, function(x) c(x$conf[,2])); poi2a <- unlist(poi2)
    yMa <- max(poi2a); yMi <- min(poi2a)
    xMa <- 1.55; xMi <- -1.25; yMa <- 1.2; yMi <- -1.05
    
    ########################
    pdf(sprintf("Output/FigS20water/Minsegs_%s/Trajectories_%s_Depth.pdf", Msegs[m], tax[i]), width=16, height=4.4)
    par(mfrow = c(1, 4))
    for(j in 1:length(W_density)){
      summ <- summary(EDR_list[[i]][[j]])
      summ <- summ[(is.na(summ$Avg_depth) == FALSE), ] #-- remove NA
      srt <- summ$ID[which.max(summ$Avg_depth)] #-- select max length
      srn <- as.numeric(str_sub(srt, 2,-1))
      
      plot(x = EDR_list[[i]][[j]], d = d_list[[i]][[j]], trajectories = t_list[[i]][[j]], states = s_list[[i]][[j]],
           select_RT = srt, #if we need to highlight some RT
           traj.colors = "lightblue", RT.colors = "orange", sel.color = "darkgreen",
           link.lty = 1, asp = 1, main = sprintf("Dilution rate : %s", W_density[j]),
           xlim = c(xMi, xMa), ylim = c(yMi, yMa), xaxs="i", yaxs="i")
    }
    dev.off()
  }
}
