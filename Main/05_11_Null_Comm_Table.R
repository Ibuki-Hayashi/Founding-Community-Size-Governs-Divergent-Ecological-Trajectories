source("Functions/functions.R")

#-- information of directory
fd <- "Output"
Sil <- readRDS(sprintf("%s/Multinomialfir_Master_ACR.rds", fd))
tax <- c("cla","ord","fam","gen","ASV","OTU99","OTU97")

df <- Sil[[6]]

#-- Preservation as excel format
if(length(list.dirs(sprintf("%s/TableS8", fd), recursive=F))==0){
  dir.create(sprintf("%s/TableS8", fd))
}

wb <- createWorkbook()
addWorksheet(wb, "Sheet1")
writeData(wb, "Sheet1", df)

saveWorkbook(wb, sprintf("%s/TableS8/Bray_Silres.xlsx", fd), overwrite = TRUE)
