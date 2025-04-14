library(vegan); library(stringr)

#-- read data from previous dir
source("Functions/functions.R")
exp <- read.csv("Data/expdata_rep.csv")
exp$Inoculum_density <- as.character(exp$Inoculum_density)
exp$Plate_position <- as.character(exp$Plate_position)

number <- 3000 #- !!!!!
Aseqtab <- readRDS(sprintf("Output/ASV_seqtab_%s.rds", as.character(number)))
Rseqtab <- readRDS(sprintf("Output/List_seqtab_%s.rds", as.character(number)))
taxtable <- readRDS(sprintf("Output/taxonomylist_%s.rds", as.character(number)))

ti <- c(2,4,6,8) # important!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! repするtimeです

#- merge_rep_data
merge_rep <- function(seqtab,exp=exp,ti=ti){
  st <- cbind(Sample_name = rownames(seqtab),as.data.frame(seqtab))
  data <- merge(exp,st,by="Sample_name"); data <- data[(data$Day %in% ti),]
  rownames(data) <- as.character(data$Sample_name)
  ds <- list(); dw <- list()
  for(i in 1:length(ti)){
    d1 <- data[(data$Day == ti[i])&(data$Source == "soil"),]
    d2 <- data[(data$Day == ti[i])&(data$Source == "water"),]
    ds[[i]] <- as.numeric(str_sub(d1$Sample_name,-4,-1))%%384
    dw[[i]] <- as.numeric(str_sub(d2$Sample_name,-4,-1))%%384
  }
  key_s <- intersect(ds[[1]],ds[[2]])
  for(i in 3:length(ti)){
    key_s <- intersect(key_s,ds[[i]])
  }
  key_w <- intersect(dw[[1]],dw[[2]])
  for(i in 3:length(ti)){
    key_w <- intersect(key_w,dw[[i]])
  }
  dd_s <- data[(data$Source == "soil"),]; dd_w <- data[(data$Source == "water"),]
  data <- rbind(dd_s[((as.numeric(str_sub(dd_s$Sample_name,-4,-1))%%384) %in% key_s),],dd_w[((as.numeric(str_sub(dd_w$Sample_name,-4,-1))%%384) %in% key_w),])
  ans <- list(data,key_s,key_w,ti); return(ans)
}
merge_d <- function(seqtab,exp=exp,ti=ti){
  st <- cbind(Sample_name = rownames(seqtab),as.data.frame(seqtab))
  data <- merge(exp,st,by="Sample_name"); data <- data[(data$Day %in% ti),]
  rownames(data) <- as.character(data$Sample_name)
  
  ans <- list(data,ti); return(ans)
}

###### return genus, family table#####
#-- set analysis dir
fd <- "Output"

out1 <- list(); out2 <- list()
for(i in 1:length(Rseqtab)){
  hoge1 <- merge_rep(seqtab=Rseqtab[[i]],exp=exp,ti=ti)
  out1[[i]] <- hoge1[[1]]
  hoge2 <- merge_d(seqtab=Rseqtab[[i]],exp=exp,ti=ti)
  out2[[i]] <- hoge2[[1]]
}

a_out <- merge_rep(seqtab=Aseqtab,exp=exp,ti=ti)
out1[[(length(Rseqtab)+1)]] <- list(a_out[[2]],a_out[[3]])
a_out2 <- merge_d(seqtab=Aseqtab,exp=exp,ti=ti)
out2[[(length(Rseqtab)+1)]] <- a_out2[[2]]

saveRDS(out1,sprintf("%s/%s_data_%s.rds",fd,"List","rep"))
saveRDS(out2,sprintf("%s/%s_data_%s.rds",fd,"List", number))
