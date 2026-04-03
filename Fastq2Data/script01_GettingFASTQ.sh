#!/bin/bash

set -e

mkdir log_and_script
mkdir -p 01_Demultiplexed_fastaq
BASE_URL="ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/fastq"

#######################################################
# The sequencing data were separately deposited in DDBJ by sequencing runs
# Downstream quality filtering was conducted on each sequencing run
# First, download_run get FASTQs included in each sequencing run (Run 1-6)
# DRA/DRX/DRR accessions are:

# [DRA accession]
# Run1: 021848
# Run2: 021850
# Run3: 021851
# Run4: 021852
# Run5: 021853
# Run6: 021854

# [DRX accession]
# Run1: 691719-692494 (the number of sample is 776)
# Run2: 692496-694031 (the number of sample is 1536)
# Run3: 694032-695567 (the number of sample is 1536)
# Run4: 695568-696335 (the number of sample is 768)
# Run5: 696336-697103 (the number of sample is 768)
# Run6: 697104-697304 (the number of sample is 201)

# [DRR accession]
# Run1: 711698-712473 (the number of sample is 776)
# Run2: 712475-714010 (the number of sample is 1536)
# Run3: 714011-715546 (the number of sample is 1536)
# Run4: 715547-716314 (the number of sample is 768)
# Run5: 716315-717082 (the number of sample is 768)
# Run6: 717083-717283 (the number of sample is 201)
#######################################################

# Downloading fastq by wget via ftp from DDBJ
download_run () {
  RUN_NAME=$1
  DRA=$2
  DRX_START=$3
  DRX_END=$4
  DRR_START=$5
  DRR_END=$6

  echo "Processing ${RUN_NAME}..."

  mkdir -p 01_Demultiplexed_fastaq/${RUN_NAME}
  cd 01_Demultiplexed_fastaq/${RUN_NAME}

  # DRA XXX (e.g. 021848 -> 021)
  DRA_PREFIX=${DRA:0:3}

  i=0
  for DRX in $(seq $DRX_START $DRX_END); do
    DRR=$(($DRR_START + i))

    DRX_ID=$(printf "DRX%06d" $DRX)
    DRR_ID=$(printf "DRR%06d" $DRR)

    URL="${BASE_URL}/DRA${DRA_PREFIX}/DRA${DRA}/${DRX_ID}/${DRR_ID}.fastq.bz2"

    echo "Downloading ${DRR_ID}..."

    wget -q $URL

    # unzip
    bunzip2 "${DRR_ID}.fastq.bz2"

    i=$((i + 1))
  done

  cd ..
}

# =========================
# for each Run
# =========================

download_run Run1 021848 691719 692494 711698 712473
download_run Run2 021850 692496 694031 712475 714010
download_run Run3 021851 694032 695567 714011 715546
download_run Run4 021852 695568 696335 715547 716314
download_run Run5 021853 696336 697103 716315 717082
download_run Run6 021854 697104 697304 717083 717283

echo "All downloads completed."


# Renaming fastq
# Each name of FASTQ is based on accession number,
# but actually they are based on the name of samples actual analysis 

cat <<'RRR' > log_and_script/script01_1_rename_fastq.R
#!/usr/bin/env Rscript

## ===================== Settings ===================== ##
csv_file <- "SraRunTable.csv"
suffix <- "515f.forward.fastq"

## ===================== Load metadata ===================== ##
df <- read.csv(csv_file, stringsAsFactors = FALSE)

if ("Sample_name" %in% colnames(df)) {
  sample_col <- "Sample_name"
} else if ("Sample Name" %in% colnames(df)) {
  sample_col <- "Sample Name"
} else {
  stop("Sample name column not found")
}

if (!"Run" %in% colnames(df)) {
  stop("Run column not found")
}

## ===================== Process ===================== ##
for (i in seq_len(nrow(df))) {

  run <- df$Run[i]
  sample <- df[[sample_col]][i]

  # File searching
  pattern <- paste0("^", run, "\\.fastq$")
  files <- list.files(".", pattern = pattern, recursive = TRUE, full.names = TRUE)

  if (length(files) == 0) {
    message("Warning: ", run, ".fastq not found")
    next
  }

  if (length(files) > 1) {
    message("Warning: multiple files found for ", run, ", using first one")
  }

  input_file <- files[1]

  output_dir <- dirname(input_file)

  output_file <- file.path(
    output_dir,
    paste0("exp_", sample, "_", suffix, ".gz")
  )

  message("Processing ", input_file, " -> ", output_file)

  # ===== gzip =====
  success <- FALSE

  tryCatch({

    con_in <- file(input_file, "rb")
    con_out <- gzfile(output_file, "wb")

    repeat {
      chunk <- readBin(con_in, what = raw(), n = 1e6)
      if (length(chunk) == 0) break
      writeBin(chunk, con_out)
    }

    close(con_in)
    close(con_out)

    success <- TRUE

  }, error = function(e) {
    message("Error during compression: ", input_file)
    message(e)
  })

  # ===== Delete original =====
  if (success) {
    file.remove(input_file)
    message("Deleted original: ", input_file)
  } else {
    message("Skipped deletion due to error: ", input_file)
  }
}

message("All done.")
RRR

Rscript script01_1_rename_fastq.R