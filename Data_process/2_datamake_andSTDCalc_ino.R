#-- information of directory
fd <- "Output"; dd <- "Data"

#-- read data from previous dir
source("Functions/functions.R")
seqtab <- readRDS(sprintf("%s/merge_ASVtab.rds", dd))
copy <- readRDS(sprintf("%s/HS3_eachASV.rds", fd))
OTU99tab <- readRDS(sprintf("%s/seqOTUtab_99.rds", fd)); class(OTU99tab) <- "numeric"
OTU97tab <- readRDS(sprintf("%s/seqOTUtab_97.rds", fd)); class(OTU97tab) <- "numeric"
taxtable <- readRDS(sprintf("%s/merge_ASVtaxonomylist.rds", dd))
sout_taxtable <- taxtable[taxtable[,3] != "c_",]
exp_i <- read.csv(sprintf("%s/expdata_ino.csv", dd))

#- soil = _5000, water = _50000 (soil's 5000 = bairitu_5000, water's 50000 = bairitu_1000)
Dsoil <- colSums(seqtab[(str_sub(rownames(seqtab), 1, 14) == "Ino_soil_5000_"),])
Dsoil_99 <- colSums(OTU99tab[(str_sub(rownames(OTU99tab), 1, 14) == "Ino_soil_5000_"),])
Dsoil_97 <- colSums(OTU97tab[(str_sub(rownames(OTU97tab), 1, 14) == "Ino_soil_5000_"),])

Dwater <- colSums(seqtab[(str_sub(rownames(seqtab), 1, 16) == "Ino_water_50000_"),])
Dwater_99 <- colSums(OTU99tab[(str_sub(rownames(OTU99tab), 1, 16) == "Ino_water_50000_"),])
Dwater_97 <- colSums(OTU97tab[(str_sub(rownames(OTU97tab), 1, 16) == "Ino_water_50000_"),])

soil <- Dsoil[(names(Dsoil) %in% rownames(taxtable[taxtable[,3] != "c_",]))]
soil_99 <- Dsoil_99[(names(Dsoil_99) %in% rownames(taxtable[taxtable[,3] != "c_",]))]
soil_97 <- Dsoil_97[(names(Dsoil_97) %in% rownames(taxtable[taxtable[,3] != "c_",]))]

water <- Dwater[(names(Dwater) %in% rownames(taxtable[taxtable[,3] != "c_",]))]
water_99 <- Dwater_99[(names(Dwater_99) %in% rownames(taxtable[taxtable[,3] != "c_",]))]
water_97 <- Dwater_97[(names(Dwater_97) %in% rownames(taxtable[taxtable[,3] != "c_",]))]

if(length(list.dirs(sprintf("%s/FigS1", fd), recursive=F))==0){
  dir.create(sprintf("%s/FigS1", fd))
}

pdf(sprintf("%s/FigS1/Inoculum_Rarefy.pdf", fd))
rarecurve(rrarefy(soil, sum(soil)), main="Soil-Inoculum", label=F)
rarecurve(rrarefy(water, sum(water)), main="Water-Inoculum", label=F)
dev.off()
####################################
#- Calculation of DNA Concentration
# -- Set random seeds (for reproduction)
ran.seed <- 1234
set.seed(ran.seed)

stds <- 5000; stdw <- 1000

STDSTD <- function(seq, taxa.print, std.bairitu, seed=1234){ #- seq=Dsoil_99; taxa.print=taxtable; std.bairitu=stds
  set.seed(seed)
  std.copy.n <- c(0.1, 0.05, 0.02, 0.01,0.005)*(6.02*10^14/1000000)/std.bairitu
  
  detected.std.name <- unique(taxa.print[which(substr(taxa.print[,"Phylum"], 1, 7) == "STD_pro"), "Phylum"])
  std.table <- t(as.data.frame(seq[names(seq) %in% rownames(taxa.print[(taxa.print[,3] == "c_"),])]))
  
  std.taxa <- c()
  for(i in 1:ncol(std.table)){
    k <- colnames(std.table)[i]
    std.taxa[i] <- taxa.print[(rownames(taxa.print) == k), 2]
  }
  
  names(std.table) <- std.taxa
  
  MergeSTD <- function(std.i, std.data = std.table){ #- std.i=detected.std.name[1]
    index.std <- which(match(names(std.table), std.i) == 1)
    if(length(index.std) > 1){
      std.tmp <- sum(std.table[index.std])
    }else{
      std.tmp <- std.table[index.std]
    }
    return(std.tmp)
  }
  new.std.table <- data.frame(std_rank1 = MergeSTD(detected.std.name[1], std.data = std.table),
                              std_rank2 = MergeSTD(detected.std.name[2], std.data = std.table),
                              std_rank3 = MergeSTD(detected.std.name[3], std.data = std.table),
                              std_rank4 = MergeSTD(detected.std.name[4], std.data = std.table),
                              std_rank5 = MergeSTD(detected.std.name[5], std.data = std.table))
  
  new.std.table[is.na(new.std.table)] <- 0
  adj.r.fun <- function(x) summary(lm(as.numeric(x) ~ std.copy.n + 0))$adj.r.squared
  lm.coef.fun <- function(x) summary(lm(as.numeric(x) ~ std.copy.n + 0))$coefficients[1]
  r2.summary <- apply(new.std.table, 1, adj.r.fun)
  coef.summary <- apply(new.std.table, 1, lm.coef.fun)
  
  new.seqtab <- seq[!(names(seq) %in% rownames(taxa.print[(taxa.print[, 3] == "c_"), ]))] # Make seq table without standard DNA
  new.seqtab2 <- new.seqtab
  
  ans <- list()
  ans[[1]] <- new.seqtab/coef.summary
  ans[[2]] <- r2.summary
  
  return(ans)
}
hosei_soil <- STDSTD(seq=Dsoil, taxa.print=taxtable, std.bairitu=stds)[[1]]
hosei_soil_99 <- STDSTD(seq=Dsoil_99, taxa.print=taxtable, std.bairitu=stds)[[1]]
hosei_soil_97 <- STDSTD(seq=Dsoil_97, taxa.print=taxtable, std.bairitu=stds)[[1]]

hosei_water <- STDSTD(seq=Dwater, taxa.print=taxtable, std.bairitu=stdw)[[1]]
hosei_water_99 <- STDSTD(seq=Dwater_99, taxa.print=taxtable, std.bairitu=stdw)[[1]]
hosei_water_97 <- STDSTD(seq=Dwater_97, taxa.print=taxtable, std.bairitu=stdw)[[1]]

###### return ans#####
#-- set analysis dir
hans_soil <- taxaratio(t(as.data.frame(hosei_soil)), sout_taxtable)
hans_water <- taxaratio(t(as.data.frame(hosei_water)), sout_taxtable)
ans_soil <- taxaratio(t(as.data.frame(soil)), sout_taxtable)
ans_water <- taxaratio(t(as.data.frame(water)), sout_taxtable)

ans_soil[[(length(ans_soil)+1)]] <- soil; ans_soil[[(length(ans_soil)+1)]] <- soil_99; ans_soil[[(length(ans_soil)+1)]] <- soil_97
ans_water[[(length(ans_water)+1)]] <- water; ans_water[[(length(ans_water)+1)]] <- water_99; ans_water[[(length(ans_water)+1)]] <- water_97
hans_soil[[(length(hans_soil)+1)]] <- hosei_soil; hans_soil[[(length(hans_soil)+1)]] <- hosei_soil_99; hans_soil[[(length(hans_soil)+1)]] <- hosei_soil_97
hans_water[[(length(hans_water)+1)]] <- hosei_water; hans_water[[(length(hans_water)+1)]] <- hosei_water_99; hans_water[[(length(hans_water)+1)]] <- hosei_water_97

#-- The number of 16S copies about each ASV(& OTU)
Conv_16S <- function(tab, Corr){ #- Corr = copy$tab
  uk <- names(tab); 
  for(i in 1:length(uk)){
    coef <- Corr[(Corr$label == uk[i]), 2]
    tab[i] <- tab[i]/coef
  }
  return(tab)
}

hs <- list(); hw <- list()
for(i in 5:7){
  hs[[i]] <- Conv_16S(tab=hans_soil[[i]], Corr=copy$tab)
  hw[[i]] <- Conv_16S(tab=hans_water[[i]], Corr=copy$tab)
}

Soil_ans <- list(ans_soil, hans_soil, hs)
Water_ans <- list(ans_water, hans_water, hw)

saveRDS(Soil_ans, sprintf("%s/Inoculum_data_soil.rds", fd))
saveRDS(Water_ans, sprintf("%s/Inoculum_data_water.rds", fd))
