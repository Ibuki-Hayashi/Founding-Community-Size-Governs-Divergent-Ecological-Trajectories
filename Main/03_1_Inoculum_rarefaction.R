source("Functions/functions.R")

#-- information of directory
fd <- "Output"; dd <- "Data"
exp <- read.csv(sprintf("%s/expdata_rep.csv", dd))
S_data <- readRDS(sprintf("%s/Inoculum_data_soil.rds", fd))
W_data <- readRDS(sprintf("%s/Inoculum_data_water.rds", fd))

sASV <- S_data[[1]][[5]]
wASV <- W_data[[1]][[5]]

# -- output
if(length(list.dirs(sprintf("%s/FigS1", fd), recursive=F))==0){
  dir.create(sprintf("%s/FigS1", fd))
}

pdf(sprintf("%s/FigS1/Inoculum_Rarefy_ASV.pdf", fd), width=8, height=4)
par(mfrow = c(1, 2))
rarecurve(x=t(sASV), step=50, label=TRUE, ylim=c(0, 150), xlim=c(0, 33000),
          col="red", main="Soil-Inoculum", ylab="The number of ASVs", xlab="Reads")
rarecurve(t(wASV), step=50, label=TRUE, ylim=c(0, 40), xlim=c(0, 4500),
          col="blue", main="Water-Inoculum", ylab="The number of ASVs", xlab="Reads")
dev.off()

