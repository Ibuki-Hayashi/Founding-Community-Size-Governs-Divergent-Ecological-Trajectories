#- Packages installation script

install.packages(c("remotes", "BiocManager", "devtools"))
library(remotes)

# --- CRAN packages (version fixed) ---
#- R version is 4.5.2

cran_pkgs <- list(
  nloptr       = "2.2.1", #- dependency of ecoregime
  lme4         = "1.1-38", #- dependency of ecoregime
  lmerTest     = "3.1-3", #- dependency of ecoregime
  jpeg         = "0.1-11", #- dependency of ggtext
  gridtext     = "0.1.5", #- dependency of ggtext
  microbiome   = "1.26.0", #- dependency of RasperGade16S
  
  cowplot      = "1.2.0",
  doParallel   = "1.0.17",
  dplyr        = "1.1.4",
  ecoregime    = "0.3.0",
  foreach      = "1.5.2",
  ggforce      = "0.5.0",
  ggnewscale   = "0.5.2",
  ggplot2      = "4.0.1",
  ggplate      = "0.2.0",
  ggsci        = "4.2.0",
  ggtext       = "0.1.2",
  ggh4x        = "0.3.1",
  gtools       = "3.9.5",
  gridExtra    = "2.3",
  igraph       = "2.2.1",
  lemon        = "0.5.2",
  lmerTest     = "3.1-3",
  lmodel2      = "1.7.4",
  markdown     = "2.0",
  MCMCpack     = "1.7-1",
  multimode    = "1.5",
  openxlsx     = "4.2.8.1",
  philentropy  = "0.10.0",
  plotly       = "4.11.0",
  plyr         = "1.8.9",
  ppcor        = "1.1",
  RColorBrewer = "1.1-3",
  reshape2     = "1.4.5",
  scales       = "1.4.0",
  seqinr       = "4.2-36",
  shortread    = "1.52.0",
  stringdist   = "0.9.15",
  smacof       = "2.1.7",
  stringr      = "1.6.0",
  tibble       = "3.3.0",
  tidyr        = "1.3.1",
  tidygraph    = "1.3.1",
  tidyverse    = "2.0.0",
  vegan        = "2.6-6.1"
)

for (pkg in names(cran_pkgs)) {
  remotes::install_version(pkg, version = cran_pkgs[[pkg]])
}

# --- Bioconductor packages ---
BiocManager::install(version = "3.22")
BiocManager::install("phyloseq")
BiocManager::install("microbiome")

# --- GitHub packages (not in CRAN) ---
devtools::install_github("wu-lab-uva/RasperGade")
devtools::install_github("wu-lab-uva/RasperGade16S")

devtools::install_github("benjjneb/dada2", ref="v1.18")

#-- rELA is not on CRAN or Bioconductor, but on GitHub
#- Please download 'rELA.v0.81.tar.gz' from https://github.com/kecosz/rELA
#- and place it in the same directory as this script, then run the following command to install it

if(F){
  install.packages("rELA/rELA.v0.81.tar.gz", type = "source")
}

