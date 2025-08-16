rm(list = ls(all.names = TRUE))

library(ggplot2); library(vegan);library(gridExtra);library(lemon);library(tidyr);library(markdown)
library(stringr);library(ggforce);library(cowplot);library(ggnewscale); library(reshape2);library(ggtext)
library(ecoregime); library(grid); library(ppcor); library(lmodel2); library(dplyr)
source("Functions/functions.R")

fd <- "Quantification_OTU"

folder <- c(list.files("Output"), "hoge")
if(any(!(folder %in% fd))){
  dir.create(sprintf("Output/%s", fd))
}
fd <- sprintf("Output/%s", fd)

dd1 <- "Data"; dd2 <- "Data2"

Daycol <- c("#1f78b4", "#33a02c", "#ff7f00", "#6a3d9a")
Scol <- c("#8B2500", "#004EbB")
Inocol <- c("#4B0082", "#800080", "#C71585", "#E6CFE6")

#-- information of directory
exp <- read.csv(sprintf("%s/expdata_rep.csv", dd1))
data <- readRDS(sprintf("%s/List_data_rep.rds", dd2)); tax <- c("cla","ord","fam","gen","ASV","OTU99","OTU97")
Sil <- readRDS(sprintf("%s/Binomialfir_Master_HY.rds", dd2))
S_density <- W_density <- c(1, 0.1, 0.01, 0.001)
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

#-- MA regression
MA_regression <- function(formula, data, weights = NULL) {
  # lmodel2(MA regression)
  model <- lmodel2(formula, data = data)
  
  # slope & intercept
  slope <- model$regression.results$Slope[2]
  intercept <- model$regression.results$Intercept[2]
  
  # returning coefficients
  return(coef(lm(y ~ x, data = data, offset = intercept + slope * data$x - data$x)))  
}


#########################################################
##-- Correlation between Continuity and Discreteness all
#########################################################

#-- Silplot
pl_A <- list(); pl_C <- list(); pl_D <- list()
Sanko <- list(); Sanko2 <- list()

for(i in 6){ #- i=6
  SilRes <- Sil[[i]]
  
  uSil <- SilRes[(SilRes$LL != Inf),]
  uSil <- SilRes[(SilRes$CV != Inf),]
  uSil$Ino_CV <- log(uSil$Ino_CV, 10)
  uSil$unique <- sprintf("%s_%s", uSil$Source, uSil$Taxa)
  uSil$Day <- as.factor(uSil$Day)
  uSil$Density <- factor(uSil$Density, levels=c(1, 0.1, 0.01, 0.001))
  uSil$CV <- (uSil$CV-min(uSil$CV))/(max(uSil$CV)-min(uSil$CV))
  uSil$Statistic <- (uSil$Statistic-min(uSil$Statistic))/(max(uSil$Statistic)-min(uSil$Statistic))
  
  pf <- function(us=uSil, dd, exp, col=Inocol, col2=Scol, leg=0){ #- dd=2; exp=1
    df <- us[us$Day == dd,]
    if(exp == 1){
      df <- data.frame(Source=df$Source, Density=df$Density,
                       Day=df$Day, x=df$Ino_CV, y=df$Statistic)
    }else{
      df <- data.frame(Source=df$Source, Density=df$Density,
                       Day=df$Day, x=df$Ino_CV, y=df$CV)
    }
    
    ma_results <- df %>%
      group_by(Source) %>%
      summarise(
        model = list(lmodel2(y ~ x, data = cur_data()))
      ) %>%
      rowwise() %>%
      mutate(
        intercept = model$regression.results$Intercept[2],
        slope = model$regression.results$Slope[2],
        
        # Calculating SE
        se_intercept=(model$confidence.intervals[2,3] - model$confidence.intervals[2,2])/(2*1.96),
        se_slope=(model$confidence.intervals[2,5] - model$confidence.intervals[2,4])/(2*1.96)
      ) 
    
    se_data <- df %>%
      group_by(Source) %>%
      summarise(
        x_min = min(x),
        x_max = max(x)
      ) %>%
      rowwise() %>%
      mutate(
        x_range = list(seq(x_min, x_max, length.out = 100)) 
      ) %>%
      unnest(x_range) %>%
      left_join(ma_results, by = "Source") %>%
      mutate(
        y_fit = intercept + slope * x_range,
        y_min = y_fit - (se_slope * x_range + se_intercept),
        y_max = y_fit + (se_slope * x_range + se_intercept)
      )
    
    line_data <- df %>%
      group_by(Source) %>%
      summarise(
        x_min = min(x),
        x_max = max(x)
      ) %>%
      left_join(ma_results, by = "Source") %>%
      mutate(
        y_min = intercept + slope * x_min,
        y_max = intercept + slope * x_max
      )
    
    ans <- ggplot(df, aes(x=x, y=y))+
      geom_ribbon(data = se_data, aes(x=x_range, y=y_fit, ymin = y_min, ymax = y_max, fill = Source), alpha = 0.16)+
      theme_bw()+scale_fill_manual(values=col2)+ggtitle(sprintf("Day%s", dd))+
      theme(text=element_text(size=7.2),
            axis.title=element_markdown(),
            legend.position = "none")+
      new_scale_fill()+
      geom_segment(data = line_data, aes(x = x_min, xend = x_max, y = y_min, yend = y_max, group = Source),
                   size = 0.3, linetype=2, color="grey50")+
      geom_point(aes(shape=Source, fill=Density), size=1.3, color="grey10", stroke=0.25)+
      scale_fill_manual(values=col)+xlab("Initial variation<br>upon inoculation")+ylim(c(-0.08, 1.1))+
      scale_shape_manual(values = c(21, 24))
    if(exp==1){
      ans <- ans+ylab("Multimodality")
    }else{
      ans <- ans+ylab("Variation caused<br>after inoculation")
    }
    
    if(leg == 1){
      ans <- ggplot(df, aes(x=x, y=y))+
        geom_ribbon(data = se_data, aes(x=x_range, y=y_fit, ymin = y_min, ymax = y_max, fill = Source), alpha = 0.16)+
        theme_bw()+scale_fill_manual(values=col2)+ggtitle(sprintf("Day%s", dd))+theme(text=element_text(size=7.2))+
        new_scale_fill()+
        geom_abline(data = ma_results, aes(slope = slope, intercept = intercept, group = Source),
                    size = 0.3, linetype=2, color="grey50")+
        geom_point(aes(shape=Source, fill=Density), size=1.3, color="grey10", stroke=0.25)+
        scale_fill_manual(values=col)+xlab("")+ylab("")+
        scale_shape_manual(values = c(21, 24))
    }
    return(ans)
  }
  c2 <- pf(us=uSil, dd=2, exp=2)
  c4 <- pf(us=uSil, dd=4, exp=2)
  c6 <- pf(us=uSil, dd=6, exp=2)
  c8 <- pf(us=uSil, dd=8, exp=2)
  
  d2 <- pf(us=uSil, dd=2, exp=1)
  d4 <- pf(us=uSil, dd=4, exp=1)
  d6 <- pf(us=uSil, dd=6, exp=1)
  d8 <- pf(us=uSil, dd=8, exp=1)
  
  cl <- pf(us=uSil, dd=2, exp=2, leg=1)
  cl <- g_legend(cl)
  
  pl_C[[i]] <- plot_grid(plot_grid(c2, c4, c6, c8, nrow=2, ncol=2), cl, ncol=2, rel_widths = c(1,0.2))
  pl_D[[i]] <- plot_grid(plot_grid(d2, d4, d6, d8, nrow=2, ncol=2), cl, ncol=2, rel_widths = c(1,0.2))
  
  pl_A[[i]] <- plot_grid(pl_D[[i]], pl_C[[i]], ncol=2)
}

pdf(sprintf("%s/ForFig_Disc_and_Cont_All2.pdf", fd), width=8, height=3.8)
pl_A[[6]]
dev.off()
