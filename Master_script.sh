#!/bin/bash
set -e

#################################################

#################################################


#- please set in the top-level directory
cd ../

###- if you need to install packages
#- Rscript Functions/install_library.R

###- Making data
echo "Step1: Making data"

Rscript Data_process/0_allcheck.R
Rscript Data_process/0_check_repcount.R
Rscript Data_process/0_copynumber.R
Rscript Data_process/0_stdout_fasta.R
Rscript Data_process/1_datamake.R
Rscript Data_process/2_datamake_andSTDCalc_ino.R
Rscript Data_process/3_datamake_rep.R

echo "Making data is completed"


###- Main scripts
echo "Step2: Main analysis"

Rscript Main/01_1_alphas_change_384_effectiveshannon.R
Rscript Main/01_2_alphas_change_384_shannonlogS.R
Rscript Main/01_3_alphas_change_384.R
Rscript Main/02_1_barplots_forFig.R
Rscript Main/02_2_barplots_supplement1.R
Rscript Main/02_3_barplots_supplement2.R
Rscript Main/03_1_Inoculum_rarefaction.R
Rscript Main/03_2_Inoculum_lott.R
Rscript Main/03_3_Inoculum_summary.R
#- Rscript Main/04_1_OTUfitting_Master.R
#- Rscript Main/04_2_OTUfitting_Master_CLR.R
#- Rscript Main/04_3_OTUfitting_Master_HH_CLR.R
Rscript Main/04_4_OTUfitting_eachOTU_plot.R
#- Rscript Main/04_5_OTUfitting_Histogram.R
#- Rscript Main/04_6_OTUfitting_Histogram_Fig2.R
#- Rscript Main/04_7_OTUfitting_Master_HHvsACR.R
Rscript Main/04_8_OTUfitting_eachOTUDay0_8.R
Rscript Main/04_9_OTUfitting_mainfig.R
Rscript Main/04_10_OTUfitting_mainfig_HH.R
Rscript Main/04_11_TableS2_GLMM.R
Rscript Main/04_12_OTUfitting_Table.R
#- Rscript Main/05_1_Null_Comm.R
#- Rscript Main/05_2_Null_Comm_HH.R
#- Rscript Main/05_3_Null_Comm_JensenShannon.R
#- Rscript Main/05_4_Null_Comm_Hellinger.R
#- Rscript Main/05_5_Null_Comm_Jaccard.R
Rscript Main/05_6_Null_Comm_Fig.R
Rscript Main/05_7_Null_Comm_Fig_HH.R
Rscript Main/05_8_Null_Comm_Fig_ACRvsHH.R
Rscript Main/05_9_NullComm_Fig_JensenShannon.R
Rscript Main/05_10_Null_Comm_Fig_Hellinger.R
Rscript Main/05_11_Null_Comm_Table.R
#- Rscript Main/05_12_Histo_Comm.R
#- Rscript Main/05_13_Histo_Comm_Fig3.R
Rscript Main/06_1_Inoculum_Extinction.R
Rscript Main/06_2_Histo_Comm_251226_DAY0vs2.R
Rscript Main/06_3_NullComm_Statistic_DAY0vs2_bray.R
Rscript Main/06_4_NullComm_Statistic_DAY0vs2_jac.R
#- Rscript Main/07_1_ELA.R 	 #- Figures are not produced but you need to save plot_ly object directly in R
#- Rscript Main/07_2_ELA_ALL.R	 #- Figures are not produced but you need to save plot_ly object directly in R
Rscript Main/08_1_Ecoregime_soil.R
Rscript Main/08_2_Ecoregime_water.R
Rscript Main/08_3_Ecoregime_metrics.R
Rscript Main/08_4_Ecoregime_soil_axis.R
Rscript Main/08_5_Ecoregime_water_axis.R
Rscript Main/08_6_Ecoregime_soil_ALL.R
Rscript Main/08_7_Ecoregime_water_ALL.R
Rscript Main/09_Wellposition.R
Rscript Main/10_Temporal_beta.R
