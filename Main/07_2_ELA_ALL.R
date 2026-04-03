set.seed(123)
n.core=6

source("Functions/functions.R")

#-- information of directory
fd <- "Output"; dd <- "Data"

exp <- read.csv(sprintf("%s/expdata_rep.csv", dd))
data <- readRDS(sprintf("%s/List_data_rep.rds", fd))
tax <- c("cla","ord","fam","gen","ASV","OTU99","OTU97")

ti <- c(2,4,6,8)
Daycol <- c("#1f78b4", "#33a02c", "#ff7f00", "#6a3d9a")
source <- c("soil", "water")
den <- c("1", "0.1", "0.01", "0.001")

d <- data[[6]] #- 99% OTU data

#-- datas
ddd <- list()
for(i in 1:2){
  ddd[[i]] <- d[(d$Source == source[i]), ]
}

# selecting taxonomic level and medium
rel <- 5/(96*4*3000) # Relative frequency threshold; if less than this value, set to 0
occ <- 5/384 # Lower limit of mean occurence; species with mean occurence lower than this value are excluded.
prun <- 1.1 # Upper limit of mean occurence; species with mean occurence higher than this value are excluded.
qth <- 10^-5

for(i in 1:2){ #- i=2; i=1 #- i=1 soil; i=2 water
  input <- ddd[[i]]
  relmat <- input[,-c(1:6)]; envmat <- input[,c(1:6)]
  eladata <- Formatting(relmat, envmat, normalize=1, c(rel, occ, prun), grouping=0, grouping_th=NULL)
  
  ocmat <- eladata[[1]]
  ocmat2 <- input[,colnames(input) %in% colnames(ocmat)]
  enmat <- eladata[[3]]
  
  lab <- paste(source[i], den[j], sep="_")

  sa <- runSA(ocmat=as.matrix(ocmat), enmat=NULL, rep=256, threads=n.core)
  gstb <- gStability(sa, ocmat, enmat=NULL, th=prun, threads=n.core)
  
  dist <- vegdist(ocmat2, method="bray")
  
  pcoa <- smacofSym(dist, ndim=2)$conf
  colnames(pcoa) <- c("mMDS1", "mMDS2")
  
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
  
  mod1 <- gam(e.realize ~ s(mMDS1, mMDS2), data= df)
  
  ##
  mds1.seq <- seq(min(df$mMDS1, na.rm=TRUE), max(df$mMDS1, na.rm=TRUE), length=70)
  mds2.seq <- seq(min(df$mMDS2, na.rm=TRUE), max(df$mMDS2, na.rm=TRUE), length=70)
  
  predfun <- function(x,y){
    newdat <- data.frame(mMDS1 = x, mMDS2=y)
    predict(mod1, newdata=newdat)
  }
  fit <- outer(mds1.seq, mds2.seq, Vectorize(predfun))
  
  ###
  ## -- Plotly
  cs <- scales::rescale(quantile(fit, probs=seq(0,1,0.25)), to=c(0,1))
  
  names(cs) <-NULL
  ans <- plot_ly(data=df, x = ~mMDS1, y= ~mMDS2, z= ~e.realize) %>% 
    add_trace(data=df, x = ~mMDS1, y= ~mMDS2, z= ~e.realize,
              type = "scatter3d", mode = "markers",
              marker = list(color = ~(envmat$Day),
                            colorscale = list(seq(0,1, length=4), 
                                              Daycol),
                            size=1.2, legendgrouptitle=list(text='e.realize', font='Arial'),
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
    
  ans #- This is a plot_ly object, so you can interact with it in RStudio's viewer pane.
  
  ans2 <- plot_ly(data=df, x = ~mMDS1, y= ~mMDS2, z= ~e.realize) %>% 
    add_trace(data=df, x = ~mMDS1, y= ~mMDS2, z= ~e.realize,
              type = "scatter3d", mode = "markers",
              marker = list(color = ~(envmat$Day),
                            colorscale = list(seq(0,1, length=4), 
                                              Daycol),
                            size=1.2, legendgrouptitle=list(text='e.realize', font='Arial'),
                            line=list(width=1,color='black')), opacity = 1)  %>% 
    add_trace(data=df, x = ~mds1.seq, y= ~mds2.seq, z= ~t(fit),
              type = "surface", showscale = F,
              hidesurface=F, opacity =0.7,
              colorscale = list(cs, c('blue','lightblue','slategray', "tan", "indianred") ),
              contours = list(z=list(show = TRUE, start = min(t(fit)), end = max(t(fit)), 
                                     usecolormap=TRUE,  size=0.7, width=3))) %>% 
    layout(scene = list(xaxis = list(nticks=10, linewidth=7, gridwidth=3, showticklabels = FALSE, title=""),
                        yaxis = list(nticks=10, linewidth=7, gridwidth =3, showticklabels = FALSE, title=""),
                        zaxis = list(nticks=10, linewidth=7, gridwidth =3, showticklabels = FALSE, title=""),
                        aspectratio = list(x = .65, y = .65, z = 0.35),
                        font='Arial'))
  ans2 #- This is a plot_ly object, so you can interact with it in RStudio's viewer pane.
}

