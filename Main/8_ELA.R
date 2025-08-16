set.seed(123)
n.core=6

library(ggplot2); library(vegan);library(gridExtra);library(lemon);library(tidyr);library(markdown)
library(stringr);library(ggforce);library(cowplot);library(ggnewscale); library(reshape2);library(ggtext)
library(ecoregime)
source("Functions/functions.R")

library("doParallel"); library('tidyverse'); library('plyr'); library('gtools')
library('ggsci'); library('igraph'); library('tidygraph'); library('RColorBrewer')
library("stringdist"); library('vegan'); library("plotly"); library("parallel")
library("dplyr"); library("mgcv"); library("rELA")
library(parallel); library(foreach); library(doParallel)

#-- information of directory
fd <- "Quantification_Comm"

folder <- c(list.files("Output"), "hoge")
if(any(!(folder %in% fd))){
  dir.create(sprintf("Output/%s", fd))
}
fd <- sprintf("Output/%s", fd)

dd1 <- "Data"; dd2 <- "Data2"

exp <- read.csv(sprintf("%s/expdata_rep.csv", dd1))
data <- readRDS(sprintf("%s/List_data_rep.rds", dd2)); tax <- c("cla","ord","fam","gen","ASV","OTU99","OTU97")

#-- read data from previous dir
ti <- c(2,4,6,8)
Daycol <- c("#1f78b4", "#33a02c", "#ff7f00", "#6a3d9a")
source <- c("soil", "water")
den <- c("1", "0.1", "0.01", "0.001")

#-- datas
ddd <- list()
for(i in 1:2){
  ddd[[i]] <- list()
  for(j in 1:4){
    d <- data[[6]]
    ddd[[i]][[j]] <- d[(d$Source == source[i])&(d$Inoculum_density == den[j]), ]
  }
}

# selecting taxonomic level and medium
rel <- 5/(96*4*3000) # Relative frequency threshold; if less than this value, set to 0
occ <- 5/384 # Lower limit of mean occurence; species with mean occurence lower than this value are excluded.
prun <- 1.1 # Upper limit of mean occurence; species with mean occurence higher than this value are excluded.
qth <- 10^-5

for(i in 1:2){ #- i=2; i=1
  for(j in 1:4){ #- j=1
    input <- ddd[[i]][[j]]
    relmat <- input[,-c(1:6)]; envmat <- input[,c(1:6)]
    eladata <- Formatting(relmat, envmat, normalize=1, c(rel, occ, prun), grouping=0, grouping_th=NULL)
    
    ocmat <- eladata[[1]]
    ocmat2 <- input[,colnames(input) %in% colnames(ocmat)]
    enmat <- eladata[[3]]
    
    lab <- paste(source[i], den[j], sep="_")

    sa <- runSA(ocmat=as.matrix(ocmat), enmat=NULL, rep=256, threads=n.core)
    gstb <- gStability(sa, ocmat, enmat=NULL, th=prun, threads=n.core)
    
    dist <- vegdist(ocmat2, method="bray")
    
    pcoa <- cmdscale(dist)
    colnames(pcoa) <- c("PCoA.1", "PCoA.2")
    
    ############################################################################
    # Preparing data
    
    df <- cbind(pcoa, gstb[[1]])
    sslist <- gstb[[4]][[4]][[2]]
    colnames(sslist)[1] <- "stable.state.id"

    ############################################################################
    # Graphic parameters
    
    n.col <- length(unique(df$stable.state.id))
    colvec <- c(brewer.pal(n = n.col, name = "Set1"))
    col1 <- colvec[1:n.col]
    
    par(oma = c(0, 0, 0, 0))
    par(mar = c(0, 0, 0, 0))
    
    #######################################################
    # Landscape plotting
    
    mod1 <- gam(e.realize ~ s(PCoA.1, PCoA.2), data= df)
    
    ##
    mds1.seq <- seq(min(df$PCoA.1, na.rm=TRUE), max(df$PCoA.1, na.rm=TRUE), length=70)
    mds2.seq <- seq(min(df$PCoA.2, na.rm=TRUE), max(df$PCoA.2, na.rm=TRUE), length=70)
    
    predfun <- function(x,y){
      newdat <- data.frame(PCoA.1 = x, PCoA.2=y)
      predict(mod1, newdata=newdat)
    }
    fit <- outer(mds1.seq, mds2.seq, Vectorize(predfun))
    
    ###
    ## -- Plotly
    cs <- scales::rescale(quantile(fit, probs=seq(0,1,0.25)), to=c(0,1))
    
    names(cs) <-NULL
    #df$color=colors
    ans <- plot_ly(data=df, x = ~PCoA.1, y= ~PCoA.2, z= ~e.realize) %>% 
      add_trace(data=df, x = ~PCoA.1, y= ~PCoA.2, z= ~e.realize,
                type = "scatter3d", mode = "markers",
                marker = list(color = ~(envmat$Day),
                              colorscale = list(seq(0,1, length=4), 
                                                Daycol),
                              size=1.6, legendgrouptitle=list(text='e.realize', font='Arial'),
                              line=list(width=1,color='black')), opacity = 1)  %>% 
      add_trace(data=df, x = ~mds1.seq, y= ~mds2.seq, z= ~t(fit),
                type = "surface", showscale = F,
                hidesurface=F, opacity =0.7,
                colorscale = list(cs, c('blue','lightblue','slategray', "tan", "indianred") ),
                contours = list(z=list(show = TRUE, start = min(t(fit)), end = max(t(fit)), 
                                       usecolormap=TRUE,  size=0.7, width=3))) %>% 
      layout(scene = list(xaxis = list(nticks=10, linewidth=7, gridwidth=3),
                           yaxis = list(nticks=10, linewidth=7, gridwidth =3),
                           zaxis = list(nticks=10, linewidth=7, gridwidth =3, tickfont = list(color = 'red')),
                           aspectratio = list(x = .65, y = .65, z = 0.35),
                           font='Arial'))
    
    ans #- this is plotly object
  }
}

