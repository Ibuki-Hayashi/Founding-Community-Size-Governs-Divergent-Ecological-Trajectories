source("Functions/functions.R")

#-- information of directory
fd <- "Output"; dd <- "Data"

#-- information of directory
exp <- read.csv(sprintf("%s/expdata_rep.csv", dd))
data <- readRDS(sprintf("%s/List_data_rep.rds", fd))
Sil <- readRDS(sprintf("%s/Binomialfir_Master_ACR_CLR.rds", fd))

tax <- c("cla","ord","fam","gen","ASV","OTU99","OTU97")
S_density <- c(1, 0.1, 0.01, 0.001)
W_density <- c(1, 0.1, 0.01, 0.001)
day <- c(2, 4, 6, 8)

#-- data
sTax <- list(); wTax <- list()
for(i in 5:7){ #- i=6
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
for(i in 6){ #- i=6
  ans <- data.frame(matrix(NA, nrow=16, ncol=11))
  colnames(ans) <- c("Source","Dataset", "Objective", "Explanetory", "SE", "SD", "tvalue", "pvalue", "R2", "R2adj", "Fvalue")
  SilRes <- Sil[[i]]
  
  uSil <- SilRes[(SilRes$LL != Inf),]
  uSil <- SilRes[(SilRes$CV != Inf),]
  uSil$Ino_CV <- log(uSil$Ino_CV, 10)
  uSil$unique <- sprintf("%s_Day%s_%s", uSil$Source, uSil$Day, uSil$Taxa)
  uSil$Density <- factor(uSil$Density, levels=c(1, 0.1, 0.01, 0.001))
  
  So <- rep(rep(c("Soil", "Water"), each=4), 2)
  Ds <- rep(c("Day2Data", "Day4Data", "Day6Data", "Day8Data"), 4)
  Ob <- c(rep("Withinmode_Variation", 8), rep("Multimodality", 8))
  Ex <- c(rep("Estimated_CV", 16))
  ans$Source <- So; ans$Dataset <- Ds; ans$Objective <- Ob; ans$Explanetory <- Ex
  
  hoge <- uSil[(uSil$Day == 8)&(uSil$Source == "Soil"),]
  hoge <- uSil[(uSil$Day == 8)&(uSil$Source == "Water"),]
  
  extf <- function(y, x){#- y=uSil$CV; x=uSil$Ino_CV
    fit <- summary(lm(y~x))
    coef <- fit$coefficients
    ans <- c(fit$coefficients[2,2], fit$coefficients[2,1], fit$coefficients[2,3], fit$coefficients[2,4],
             fit$r.squared, fit$adj.r.squared, fit$fstatistic[1])
    return(ans)
  }
  
  for(j in 1:nrow(ans)){ #- j=1
    uSil2 <- uSil[uSil$Day == as.numeric(str_sub(ans[j,2], 4,4)),]
    if(ans[j,1] == "Soil"){
      uSil2 <- uSil2[uSil2$Source == "Soil",]
    }else{
      uSil2 <- uSil2[uSil2$Source == "Water",]
    }
    if(ans[j,3] == "Withinmode_Variation"){
      df <- data.frame(y=uSil2$CV, x=uSil2$Ino_CV)
    }else{
      df <- data.frame(y=uSil2$Statistic, x=uSil2$Ino_CV)
    }
    ans[j,5:11] <- extf(df$y, df$x)
  }
  
  ans <- data.frame(ans, Qvalue=p.adjust(as.numeric(ans$pvalue), method="BH"))
  
  xls[[i-5]] <- ans
}
names(xls) <- tax[6]

#-- Preservation as excel format
if(length(list.dirs(sprintf("%s/TableS7", fd), recursive=F))==0){
  dir.create(sprintf("%s/TableS7", fd))
}

wb <- createWorkbook()
for (sheet_name in names(xls)) {
  addWorksheet(wb, sheet_name)  # add sheet
  writeData(wb, sheet_name, xls[[sheet_name]])  # write data
}
saveWorkbook(wb, sprintf("%s/TableS7/MetricsVS_ACR.xlsx", fd), overwrite = TRUE)

#-- Silplot
xls <- list()
for(i in 6){ #- i=6
  ans <- data.frame(matrix(NA, nrow=12, ncol=6))
  colnames(ans) <- c("Dataset1", "Dataset2", "Category", "pvalue", "qvalue", "Significance")
  SilRes <- Sil[[i]]
  
  uSil <- SilRes[(SilRes$LL != Inf),]
  uSil <- SilRes[(SilRes$CV != Inf),]
  uSil$Ino_CV <- log(uSil$Ino_CV, 10)
  uSil$unique <- sprintf("%s_Day%s_%s", uSil$Source, uSil$Day, uSil$Taxa)
  uSil$Density <- factor(uSil$Density, levels=c(1, 0.1, 0.01, 0.001))

  xls[[i-5]] <- uSil
}
names(xls) <- tax[6]

#-- Preservation as excel format
if(length(list.dirs(sprintf("%s/TableS4", fd), recursive=F))==0){
  dir.create(sprintf("%s/TableS4", fd))
}

wb <- createWorkbook()
for (sheet_name in names(xls)) {
  addWorksheet(wb, sheet_name)  # add sheet
  writeData(wb, sheet_name, xls[[sheet_name]])  # write data
}
saveWorkbook(wb, sprintf("%s/TableS4/OTUstoc_ACR_CLR.xlsx", fd), overwrite = TRUE)
