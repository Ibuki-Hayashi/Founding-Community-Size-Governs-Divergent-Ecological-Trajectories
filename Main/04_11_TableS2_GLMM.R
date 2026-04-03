source("Functions/functions.R")

#-- information of directory
fd <- "Output"; dd <- "Data"

#-- information of directory
exp <- read.csv(sprintf("%s/expdata_rep.csv", dd))
data <- readRDS(sprintf("%s/List_data_rep.rds", fd))
Sil <- readRDS(sprintf("%s/Binomialfir_Master_ACR.rds", fd))
tax <- c("cla","ord","fam","gen","ASV","OTU99","OTU97")

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
i <- 6
ans <- data.frame(matrix(NA, nrow=16, ncol=11))
colnames(ans) <- c("Source","Dataset", "Objective", "Explanetory", "SE", "SD", "tvalue", "pvalue", "R2", "R2adj", "Fvalue")
SilRes <- Sil[[i]]
uSil <- SilRes[(SilRes$LL != Inf),]
uSil <- SilRes[(SilRes$CV != Inf),]
uSil$unique <- sprintf("%s_Day%s_%s", uSil$Source, uSil$Day, uSil$Taxa)
uSil$Density <- factor(uSil$Density, levels=c(1, 0.1, 0.01, 0.001))
  
uSil$Density <- ordered(uSil$Density, levels=c(1, 0.1, 0.01, 0.001))
uSil$Taxa <- as.factor(uSil$Taxa)
  
uSils <- uSil[(uSil$Source == "Soil")&(uSil$Day == 2),]
uSilw <- uSil[(uSil$Source == "Water")&(uSil$Day == 2),]
  
hoges <- lmerTest::lmer(uSils$Ino_CV ~ uSils$Density+(1|uSils$Taxa))
a <- summary(hoges)
coef_dfa <- as.data.frame(a$coefficients)
coef_dfa$term <- rownames(coef_dfa)

hogew <- lmerTest::lmer(uSilw$Ino_CV ~ uSilw$Density+(1|uSilw$Taxa))
b <- summary(hogew)
coef_dfb <- as.data.frame(b$coefficients)
coef_dfb$term <- rownames(coef_dfb)

xls <- list(coef_dfa, coef_dfb)
names(xls) <- c("a1", "b1")

#-- Preservation as excel format
if(length(list.dirs(sprintf("%s/TableS2", fd), recursive=F))==0){
  dir.create(sprintf("%s/TableS2", fd))
}

wb <- createWorkbook()
for (sheet_name in names(xls)) { #- sheet_name="a1"
  addWorksheet(wb, sheet_name)  # adding sheet
  writeData(wb, sheet_name, xls[[sheet_name]])  # writing data
}
saveWorkbook(wb, sprintf("%s/TableS2/GLMM.xlsx", fd), overwrite = TRUE)
