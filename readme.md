# Codes for the paper titled: Stochastic Forces in Microbial Community Assembly: Founding Community Size Governs Divergent Ecological Trajectories

This repository contains the code and processed data used in the study:

**_Stochastic Forces in Microbial Community Assembly: Founding Community Size Governs Divergent Ecological Trajectories_**  
https://doi.org/10.1101/2025.08.09.669462

---

## Overview

The repository includes:

- Processed sequencing data (post–DNA sequence preprocessing)
- Scripts for preprocessing (FASTQ to ASV/OTU tables)
- Scripts for downstream analyses

Raw sequencing data are not included. They can be obtained from:

- **DDBJ BioProject:** PRJDB35809

A download pipeline is provided but commented out due to long runtime.

---

## Workflow

The analysis consists of two main steps:

1. Preprocessing (FASTQ → ASV/OTU tables)
2. Main analysis (R-based)

---

## 1. Preprocessing

### Environment

See `Preprocess_Env.txt` for required dependencies (Unix tools and R).

### Execution

```bash
cd Fastq2Data
bash FASTQscript.sh
```

### Notes

- Data download and quality filtering are time-consuming.
- Users are encouraged to skip this step and use the provided `Data/` directory.

### Reference database

Taxonomic annotation requires:

- SILVA nr99 v138 database

Before use:

1. Append `Fastq2Data/Ref/STD.fa` to the beginning of the database file
2. Place the modified file in `Fastq2Data/Ref/`

### Output

- `Data/`

---

## 2. Main Analysis

### Environment

- Package installation scripts: `Functions/*.R`
- Session information: `sessioninfo.txt`

### Execution

```bash
bash MASTER_script.sh
```

### Notes

- Computationally intensive steps are partially commented out
- Precomputed results are included for reproducibility

### Output

- `Output/`

---

## Recommended Usage

For most users, start from the main analysis:

```bash
bash MASTER_script.sh
```

---

## Reproducibility

To reproduce the full workflow from raw data:

1. Install dependencies (`Preprocess_Env.txt`, `sessioninfo.txt`)
2. Run preprocessing
3. Run main analysis

---
