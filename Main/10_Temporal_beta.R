source("Functions/functions.R")

#-- information of directory
fd <- "Output"; dd1 <- "Data"

exp <- read.csv(sprintf("%s/expdata_rep.csv", dd))
data <- readRDS(sprintf("%s/List_data_rep.rds", fd))
tax <- c("cla","ord","fam","gen","ASV","OTU99","OTU97")
Density <- c(1, 0.1, 0.01, 0.001)
ti <- c(2, 4, 6, 8)
Daycol <- c("#1f78b4", "#33a02c", "#ff7f00", "#6a3d9a")

# -- Caluculate beta
Make_df <- function(data){ #- data = data[[6]]
  datas <- data[data$Source == "soil",]; datas$Inoculum_density <- factor(datas$Inoculum_density, levels=c("1","0.1","0.01","0.001"))
  dataw <- data[data$Source == "water",]; dataw$Inoculum_density <- factor(dataw$Inoculum_density, levels=c("1","0.1","0.01","0.001"))
  
  Ans_soil1 <- Ans_soil2 <- Ans_soil3 <- data.frame(Location = datas$Location[datas$Day == 2],
                    Inoculum_density = datas$Inoculum_density[datas$Day == 2],
                    Which = rep(NA, length(unique(datas$Location))),
                    Bray = rep(NA, length(unique(datas$Location))),
                    Jacc = rep(NA, length(unique(datas$Location))))
  Ans_water1 <- Ans_water2 <- Ans_water3 <- data.frame(Location = dataw$Location[dataw$Day == 2],
                         Inoculum_density = dataw$Inoculum_density[dataw$Day == 2],
                         Which = rep(NA, length(unique(dataw$Location))),
                         Bray = rep(NA, length(unique(dataw$Location))),
                         Jacc = rep(NA, length(unique(dataw$Location))))
  
  Ans_soil1$Which <- "Day2vs4"; Ans_soil2$Which <- "Day4vs6"; Ans_soil3$Which <- "Day6vs8"
  Ans_water1$Which <- "Day2vs4"; Ans_water2$Which <- "Day4vs6"; Ans_water3$Which <- "Day6vs8"
  
  for(i in 1:nrow(Ans_water1)){ #- i=1
    data24 <- rbind(x1=datas[(datas$Location == Ans_water1[i,1])&(datas$Day == 2), c(7:ncol(datas))],
                    x2=datas[(datas$Location == Ans_water1[i,1])&(datas$Day == 4), c(7:ncol(datas))])
    data46 <- rbind(x1=datas[(datas$Location == Ans_water2[i,1])&(datas$Day == 4), c(7:ncol(datas))],
                    x2=datas[(datas$Location == Ans_water2[i,1])&(datas$Day == 6), c(7:ncol(datas))])
    data68 <- rbind(x1=datas[(datas$Location == Ans_water3[i,1])&(datas$Day == 6), c(7:ncol(datas))],
                    x2=datas[(datas$Location == Ans_water3[i,1])&(datas$Day == 8), c(7:ncol(datas))])
    Ans_soil1[i,4] <- as.numeric(vegdist(data24)); Ans_soil1[i,5] <- as.numeric(vegdist(forjac(data24), method="jaccard"))
    Ans_soil2[i,4] <- as.numeric(vegdist(data46)); Ans_soil2[i,5] <- as.numeric(vegdist(forjac(data46), method="jaccard"))
    Ans_soil3[i,4] <- as.numeric(vegdist(data68)); Ans_soil3[i,5] <- as.numeric(vegdist(forjac(data68), method="jaccard"))
    
    data24 <- rbind(x1=dataw[(dataw$Location == Ans_water1[i,1])&(dataw$Day == 2), c(7:ncol(dataw))],
                    x2=dataw[(dataw$Location == Ans_water1[i,1])&(dataw$Day == 4), c(7:ncol(dataw))])
    data46 <- rbind(x1=dataw[(dataw$Location == Ans_water2[i,1])&(dataw$Day == 4), c(7:ncol(dataw))],
                    x2=dataw[(dataw$Location == Ans_water2[i,1])&(dataw$Day == 6), c(7:ncol(dataw))])
    data68 <- rbind(x1=dataw[(dataw$Location == Ans_water3[i,1])&(dataw$Day == 6), c(7:ncol(dataw))],
                    x2=dataw[(dataw$Location == Ans_water3[i,1])&(dataw$Day == 8), c(7:ncol(dataw))])
    Ans_water1[i,4] <- as.numeric(vegdist(data24)); Ans_water1[i,5] <- as.numeric(vegdist(forjac(data24), method="jaccard"))
    Ans_water2[i,4] <- as.numeric(vegdist(data46)); Ans_water2[i,5] <- as.numeric(vegdist(forjac(data46), method="jaccard"))
    Ans_water3[i,4] <- as.numeric(vegdist(data68)); Ans_water3[i,5] <- as.numeric(vegdist(forjac(data68), method="jaccard"))
  }
  
  Ans_soil <- rbind(Ans_soil1, rbind(Ans_soil2, Ans_soil3))
  Ans_water <- rbind(Ans_water1, rbind(Ans_water2, Ans_water3))
  
  return(list(Ans_soil, Ans_water))
}

#- make boxplot
Boxfunc <- function(df, title=""){ #- df = Make_df(data[[6]])[[1]]
  bray_p <- ggplot(df)+
    geom_boxplot(aes(x=Inoculum_density, y=Bray, fill=Which))+theme_bw()+ggtitle(title)+
    ylim(c(0,1))+scale_fill_manual(values=Daycol[1:3])+ylab("Bray-Curtis dissimilarity")
  jacc_p <- ggplot(df)+
    geom_boxplot(aes(x=Inoculum_density, y=Jacc, fill=Which))+theme_bw()+ggtitle(title)+
    ylim(c(0,1))+scale_fill_manual(values=Daycol[1:3])+ylab("Jaccard dissimilarity")
  
  ans <- plot_grid(bray_p, jacc_p, nrow=2)
  return(ans)
}

Ans_plot <- list()
for(i in 6){ #- i=6
  ss <- Boxfunc(df=Make_df(data[[i]])[[1]], title=sprintf("%s_%s", tax[i], "Soil"))
  ww <- Boxfunc(df=Make_df(data[[i]])[[2]], title=sprintf("%s_%s", tax[i], "Water"))
  
  Ans_plot[[i]] <- plot_grid(ss, ww, nrow=1)
}

# -- output
if(length(list.dirs(sprintf("%s/FigS26", fd), recursive=F))==0){
  dir.create(sprintf("%s/FigS26", fd))
}

pdf(sprintf("%s/FigS26/Temporal_beta.pdf", fd),width=7.2,height=5) #- OTU99% only
Ans_plot[[6]]
dev.off()
