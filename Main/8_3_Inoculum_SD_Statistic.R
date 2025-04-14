rm(list = ls(all.names = TRUE))

library(ggplot2); library(vegan);library(gridExtra);library(lemon);library(tidyr);library(markdown)
library(stringr);library(ggforce);library(cowplot);library(ggnewscale); library(reshape2);library(ggtext)
library(ecoregime); library(grid); library(ppcor); library(openxlsx)
source("Functions/functions.R")

fd <- "Output"; dd <- "Data"

#-- information of directory
exp <- read.csv(sprintf("%s/expdata_rep.csv", dd))
data <- readRDS(sprintf("%s/List_data_rep.rds", fd)); tax <- c("cla","ord","fam","gen","ASV","OTU99","OTU97")
Sil <- readRDS(sprintf("%s/Binomialfir_Master_HY.rds", fd))
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
xls <- list()
for(i in 5:7){ #- i=6
  ans <- data.frame(matrix(NA, nrow=10, ncol=10))
  colnames(ans) <- c("Dataset", "Objective", "Explanetory", "SE", "SD", "tvalue", "pvalue", "R2", "R2adj", "Fvalue")
  SilRes <- Sil[[i]]
  
  uSil <- SilRes[(SilRes$LL != Inf),]
  uSil <- SilRes[(SilRes$CV != Inf),]
  uSil$Ino_CV <- log(uSil$Ino_CV, 10)
  uSil$unique <- sprintf("%s_Day%s_%s", uSil$Source, uSil$Day, uSil$Taxa)
  uSil$Density <- factor(uSil$Density, levels=c(1, 0.1, 0.01, 0.001))
  
  Ds <- rep(c("FullData", "Day2Data", "Day4Data", "Day6Data", "Day8Data"), 2)
  Ob <- c(rep("Continuity", 5), rep("Discreteness", 5))
  Ex <- c(rep("CV_Inoculum", 10))
  ans$Dataset <- Ds; ans$Objective <- Ob; ans$Explanetory <- Ex
  
  extf <- function(y, x){#- y=uSil$CV; x=uSil$Ino_CV
    fit <- summary(lm(y~x))
    coef <- fit$coefficients
    ans <- c(fit$coefficients[2,2], fit$coefficients[2,1], fit$coefficients[2,3], fit$coefficients[2,4],
             fit$r.squared, fit$adj.r.squared, fit$fstatistic[1])
    return(ans)
  }
  
  for(j in 1:nrow(ans)){
    if(ans[j,1] == "FullData"){
      uSil2 <- uSil
    }else{
      uSil2 <- uSil[uSil$Day == as.numeric(str_sub(ans[j,1], 4,4)),]
    }
    if(ans[j,2] == "Continuity"){
      df <- data.frame(y=uSil2$CV, x=uSil2$Ino_CV)
    }else{
      df <- data.frame(y=uSil2$Statistic, x=uSil2$Ino_CV)
    }
    ans[j,4:10] <- extf(df$y, df$x)
  }
  xls[[i-4]] <- ans
}
names(xls) <- tax[5:7]

#-- Preservation as excel format
wb <- createWorkbook()
for (sheet_name in names(xls)) {
  addWorksheet(wb, sheet_name)  # adding a sheet
  writeData(wb, sheet_name, xls[[sheet_name]])  # writing data
}
saveWorkbook(wb, sprintf("%s/ContDist_statistic1.xlsx", fd), overwrite = TRUE)

#-- Silplot
xls <- list()
for(i in 5:7){ #- i=6
  ans <- data.frame(matrix(NA, nrow=12, ncol=6))
  colnames(ans) <- c("Dataset1", "Dataset2", "Category", "pvalue", "qvalue", "Significance")
  SilRes <- Sil[[i]]
  
  uSil <- SilRes[(SilRes$LL != Inf),]
  uSil <- SilRes[(SilRes$CV != Inf),]
  uSil$Ino_CV <- log(uSil$Ino_CV, 10)
  uSil$unique <- sprintf("%s_Day%s_%s", uSil$Source, uSil$Day, uSil$Taxa)
  uSil$Density <- factor(uSil$Density, levels=c(1, 0.1, 0.01, 0.001))
  
  Ds <- c("Day2Data", "Day4Data", "Day6Data", "Day8Data")
  Us <- t(combn(Ds, 2))
  Ob <- c(rep("Continuity", nrow(Us)), rep("Discreteness", nrow(Us)))
  ans[,c(1,2)] <-rbind(Us, Us); ans$Category <- Ob
  
  for(j in 1:nrow(ans)){ #- j=1
    if(ans[j,3] == "Continuity"){
      df1 <- uSil[(uSil$Day == as.numeric(str_sub(ans[j,1], 4,4))),]$CV
      df2 <- uSil[(uSil$Day == as.numeric(str_sub(ans[j,2], 4,4))),]$CV
      ans[j,4] <- wilcox.test(df1, df2)$p.value
    }else{
      df1 <- uSil[(uSil$Day == as.numeric(str_sub(ans[j,1], 4,4))),]$Statistic
      df2 <- uSil[(uSil$Day == as.numeric(str_sub(ans[j,2], 4,4))),]$Statistic
      ans[j,4] <- wilcox.test(df1, df2)$p.value
    }
  }
  ans[,5] <- p.adjust(ans[,4], method="BH")
  ans[,6] <- ifelse(ans[,5] < 0.01, "Sig.", "N.S.")
  xls[[i-4]] <- ans
}
names(xls) <- tax[5:7]

#-- Preservation as excel format
wb <- createWorkbook()
for (sheet_name in names(xls)) {
  addWorksheet(wb, sheet_name)  # adding a sheet
  writeData(wb, sheet_name, xls[[sheet_name]])  # writing data
}
saveWorkbook(wb, sprintf("%s/ContDist_statistic2.xlsx", fd), overwrite = TRUE)