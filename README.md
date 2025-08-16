This repository contains the code and raw data (after DNA sequence processing) used for the analysis in the paper titled *"Stochastic Forces in Microbial Community Assembly: Founding Community Size Governs Divergent Ecological Trajectories"* (bioRxiv version: https://doi.org/10.1101/2025.08.09.669462).  
The raw DNA sequences can be downloaded from DDBJ (BioProject: PRJDB35809).

1. For each of the six DRA accession numbers, the downloaded sequences are processed independently with **Script1–3** in the folder **Fastq2nochim**, which produces a file named `stall_no_rmchimera.rds` (already included in the folder **DRA**) for each DRA accession.

2. By gathering these files into a single folder (**DRA**) and running `mergeSequenceTable.R` on that folder, you can obtain the files in the **Table** folder.

3. The folder **Data_process** contains scripts that generate the files in the **Data** folder, based on the files in the **Table** folder.

4. Finally, please run the scripts in **Main** sequentially, using the files in the **Data** folder as input.
