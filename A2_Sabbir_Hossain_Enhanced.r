# ---
# A2_Sabbir_Hossain_Enhanced.R
#
# Enhanced R Script for Assignment 2: Differential Gene Expression and Preliminary ORA
# BCB420 - Computational Systems Biology
# Original Author: Sabbir Hossain
# Enhanced Version: AI Assistant
# Date: Sys.Date()
# ---

# ---
# SECTION 0: PREAMBLE AND PACKAGE MANAGEMENT
# ---
# This section ensures all necessary packages are installed and loaded.
# It's designed to be run once or checked at the beginning of each session.
# ---

cat("--- STARTING SECTION 0: PREAMBLE AND PACKAGE MANAGEMENT ---\n")

# Package installation (if needed)
# Suppressing warnings for cleaner output during package checks/installations
suppressWarnings({
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }

  # List of required Bioconductor packages
  bioc_packages <- c(
    "GEOmetadb", "limma", "AnnotationDbi", "org.Hs.eg.db",
    "edgeR", "Biobase", "biomaRt", "GEOquery"
  )

  # List of required CRAN packages
  cran_packages <- c(
    "readxl", "RColorBrewer", "tidyverse", "ggplot2",
    "dplyr", "vegan", "gprofiler2", "ComplexHeatmap", "ggrepel",
    "magrittr", "RSQLite"
  ) # Added magrittr and RSQLite explicitly

  # Install Bioconductor packages if not already installed
  for (pkg in bioc_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      BiocManager::install(pkg, update = FALSE, ask = FALSE)
    }
  }

  # Install CRAN packages if not already installed
  for (pkg in cran_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg)
    }
  }
})

# Load packages, suppressing startup messages for a cleaner console
cat("Loading required packages...\n")
suppressPackageStartupMessages({
  library(tidyverse)
  library(edgeR)
  library(GEOmetadb)
  library(RColorBrewer)
  library(ggplot2)
  library(readxl)
  library(dplyr)
  library(AnnotationDbi)
  library(limma)
  library(Biobase)
  # library(BiocManager) # Not typically loaded directly after installation
  library(biomaRt)
  library(magrittr)
  library(GEOquery)
  library(RSQLite)
  # library(org.Hs.eg.db) # Often used implicitly by AnnotationDbi or explicitly when needed
  library(vegan)
  library(gprofiler2)
  library(ComplexHeatmap) # Explicitly load for Heatmap function
  library(ggrepel)
})

cat("All required packages checked and loaded.\n")
cat("--- COMPLETED SECTION 0 ---\n\n")

# ---
# PROJECT OVERVIEW (from README.md and Rmd)
# ---
# Objective: Analyze normalized gene expression data, rank genes by differential
# expression, and conduct thresholded over-representation analysis (ORA) to
# identify key biological themes. This builds upon Assignment 1.
#
# Original Dataset: GSE193417
#   - Study on corticotropin-releasing hormone (CRH) mRNA in MDD patients.
#   - RNA-sequencing on sgACC CRH+ interneurons.
#   - Comparison between MDD participants (n=6) and controls (n=6).
# ---

# ---
# SECTION 1: DATA SETUP AND PREPROCESSING (Adapted from Assignment 1 logic)
# ---
# This section downloads, cleans, filters, and normalizes the gene expression data.
# Intermediate RDS files are used to save progress and allow faster re-runs.
# ---
cat("--- STARTING SECTION 1: DATA SETUP AND PREPROCESSING ---\n")

# Create a directory for data if it doesn't exist
if (!dir.exists("A2_data_output")) {
  dir.create("A2_data_output")
}
data_dir <- "A2_data_output/"

# File paths for RDS objects
rds_GSE193417_GEO <- file.path(data_dir, "GSE193417_GEOobject.rds")
rds_GSE193417_filtered_CPM <- file.path(data_dir, "GSE193417_filtered_CPM.rds")
rds_GSE193417_filtered_symbols <- file.path(data_dir, "GSE193417_filtered_symbols.rds")
rds_GSE193417_finalfiltered_symbols <- file.path(data_dir, "GSE193417_finalfiltered_symbols.rds")
rds_GSE193417_normalized_CPM <- file.path(data_dir, "GSE193417_normalized_CPM.rds")
rds_GSE193417_normalized_datastruct <- file.path(data_dir, "GSE193417_normalized_datastruct.rds")
rds_geneIDs <- file.path(data_dir, "GSE193417_geneIDs.rds")


# 1.1: Download GEO data (GSE193417)
cat("Step 1.1: Downloading or loading GSE193417 GEO object...\n")
if (!file.exists(rds_GSE193417_GEO)) {
  cat("  Downloading GSE193417 from GEO...\n")
  GSE193417_geo_obj_list <- getGEO("GSE193417", GSEMatrix = TRUE, getGPL = FALSE)
  if (length(GSE193417_geo_obj_list) > 1) {
    idx <- grep("GPL16791", attr(GSE193417_geo_obj_list, "names")) # Platform for this study
  } else {
    idx <- 1
  }
  GSE193417s <- GSE193417_geo_obj_list[[idx]] # This is the ExpressionSet
  saveRDS(GSE193417s, rds_GSE193417_GEO) # Saving the ExpressionSet
  cat("  GSE193417 ExpressionSet downloaded and saved.\n")
} else {
  cat("  Loading GSE193417 ExpressionSet from local RDS file...\n")
  GSE193417s <- readRDS(rds_GSE193417_GEO)
}
cat("  GSE193417 ExpressionSet object ('GSE193417s') details:\n")
print(GSE193417s)
cat("  Phenotypic data (pData) head:\n")
print(head(pData(GSE193417s)))


# 1.2: Download supplementary files (raw counts)
cat("\nStep 1.2: Downloading or locating supplementary files (raw counts)...\n")
gse_supp_files_dir <- file.path(data_dir, "GSE193417_supp")
if (!dir.exists(gse_supp_files_dir)) {
  dir.create(gse_supp_files_dir)
  cat("  Downloading supplementary files to:", gse_supp_files_dir, "\n")
  # It seems the original RMD used a local path construction after getGEOSuppFiles
  # We will download into the dedicated directory.
  getGEOSuppFiles("GSE193417", baseDir = data_dir, makeDirectory = TRUE, fetch_files = TRUE)
  # The actual file name is "GSE193417_raw_counts_GRCh38.p13_NCBI.tsv.gz"
  # getGEOSuppFiles creates a GSE193417 folder inside baseDir.
  # The path will be data_dir/GSE193417/GSE193417_raw_counts_GRCh38.p13_NCBI.tsv.gz
} else {
  cat("  Supplementary files directory", gse_supp_files_dir, "already exists. Assuming files are present.\n")
}

# Define the path to the raw counts file
# The file is GSE193417_raw_counts_GRCh38.p13_NCBI.tsv.gz
# It will be in data_dir/GSE193417/
raw_counts_file_path <- file.path(gse_supp_files_dir, "GSE193417_raw_counts_GRCh38.p13_NCBI.tsv.gz")

if (!file.exists(raw_counts_file_path)) {
  cat("  Raw counts file not found at:", raw_counts_file_path, "\n")
  cat("  Attempting to list files downloaded by getGEOSuppFiles to find it...\n")
  # This part is tricky if the exact file name isn't known or if getGEOSuppFiles behavior changes.
  # For GSE193417, there's one main count file.
  # Let's assume it's the first .tsv.gz file if an exact match fails.
  possible_files <- list.files(gse_supp_files_dir, pattern = "\\.tsv\\.gz$", full.names = TRUE)
  if (length(possible_files) > 0) {
    raw_counts_file_path <- possible_files[1]
    cat("  Found potential raw counts file:", raw_counts_file_path, "\n")
  } else {
    stop("Raw counts file could not be located. Please check download or path.")
  }
}

# 1.3: Read and initially filter raw counts
cat("\nStep 1.3: Reading and initially filtering raw counts...\n")
if (!file.exists(rds_GSE193417_filtered_CPM)) {
  cat("  Reading raw counts from:", raw_counts_file_path, "\n")
  # The RMD used read.csv, which works for .gz if recognized. read_tsv is often better for tsv.
  mddInMe <- readr::read_tsv(raw_counts_file_path, col_types = readr::cols(.default = "d", Geneid = "c"))

  # The first column in the provided file is 'Geneid' (Ensembl IDs), not EntrezID as in Rmd.
  # Let's rename it to match the Rmd's expectation for 'EntrezID' even if it's Ensembl.
  colnames(mddInMe)[1] <- "Ensembl_ID"
  cat("  Raw counts data ('mddInMe') dimensions (pre-filtering):", dim(mddInMe)[1], "genes,", dim(mddInMe)[2], "samples+ID_col.\n")
  print(head(mddInMe[, 1:4]))

  # Filter genes: keep rows where at least 6 samples have counts > 1
  # The Rmd logic: rowSums(mddInMe > 1) >= 6. This assumes mddInMe column 1 is ID.
  # The actual sample columns start from column 2.
  numeric_counts_matrix <- as.matrix(mddInMe[, -1]) # Exclude the Geneid column
  keep_genes_filter <- rowSums(numeric_counts_matrix > 1) >= 6

  filtMddInMe <- mddInMe[keep_genes_filter, ]
  cat("  Filtered counts data ('filtMddInMe') dimensions:", dim(filtMddInMe)[1], "genes,", dim(filtMddInMe)[2], "samples+ID_col.\n")
  print(head(filtMddInMe[, 1:4]))

  saveRDS(filtMddInMe, rds_GSE193417_filtered_CPM)
  cat("  Filtered counts (filtMddInMe) saved to RDS.\n")
} else {
  cat("  Loading filtered counts (filtMddInMe) from local RDS file...\n")
  filtMddInMe <- readRDS(rds_GSE193417_filtered_CPM)
  cat("  Loaded 'filtMddInMe' dimensions:", dim(filtMddInMe)[1], "genes,", dim(filtMddInMe)[2], "samples+ID_col.\n")
}

# 1.4: Map Ensembl IDs to HGNC symbols using biomaRt
cat("\nStep 1.4: Mapping Ensembl IDs to HGNC symbols...\n")
if (!file.exists(rds_geneIDs)) {
  cat("  Setting up biomaRt to connect to Ensembl...\n")
  ensembl <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  cat("  Querying Ensembl for HGNC symbols for", nrow(filtMddInMe), "filtered Ensembl IDs...\n")

  geneIDs <- getBM(
    attributes = c("ensembl_gene_id", "hgnc_symbol"),
    filters = "ensembl_gene_id",
    values = filtMddInMe$Ensembl_ID, # Use the correct column name
    mart = ensembl
  )

  cat("  biomaRt query returned", nrow(geneIDs), "mappings.\n")
  cat("  Number of unique Ensembl IDs in query:", length(unique(filtMddInMe$Ensembl_ID)), "\n")
  cat("  Number of unique Ensembl IDs in mapping result:", length(unique(geneIDs$ensembl_gene_id)), "\n")
  cat("  Number of unique HGNC symbols in mapping result:", length(unique(geneIDs$hgnc_symbol)), "\n")
  print(head(geneIDs))
  saveRDS(geneIDs, rds_geneIDs)
  cat("  Gene ID mappings (geneIDs) saved to RDS.\n")
} else {
  cat("  Loading gene ID mappings (geneIDs) from local RDS file...\n")
  geneIDs <- readRDS(rds_geneIDs)
  cat("  Loaded 'geneIDs' with", nrow(geneIDs), "mappings.\n")
}

# 1.5: Further filtering and data cleaning based on mapping
cat("\nStep 1.5: Further filtering and data cleaning based on mappings...\n")
if (!file.exists(rds_GSE193417_finalfiltered_symbols)) {
  cat("  Performing inner join between filtered counts and gene symbols...\n")
  # Ensure column names for joining are correct
  FinalGeneFilter <- dplyr::inner_join(geneIDs, filtMddInMe, by = c("ensembl_gene_id" = "Ensembl_ID"))
  cat("  Dimensions after inner join ('FinalGeneFilter'):", dim(FinalGeneFilter)[1], "genes,", dim(FinalGeneFilter)[2], "columns.\n")

  # Remove genes with empty HGNC symbols
  num_empty_hgnc_before <- sum(FinalGeneFilter$hgnc_symbol == "" | is.na(FinalGeneFilter$hgnc_symbol))
  FinalGeneFilter <- FinalGeneFilter[!(FinalGeneFilter$hgnc_symbol == "" | is.na(FinalGeneFilter$hgnc_symbol)), ]
  cat("  Removed", num_empty_hgnc_before - sum(FinalGeneFilter$hgnc_symbol == "" | is.na(FinalGeneFilter$hgnc_symbol)), "genes with empty HGNC symbols.\n")
  cat("  Dimensions after removing empty HGNC:", dim(FinalGeneFilter)[1], "genes.\n")

  # Handle duplicated HGNC symbols (keep first mapping by Ensembl ID, then by HGNC if still duplicated)
  # Or, a simpler RMD approach: remove specific problematic symbols if known
  # The RMD removed STRA6LP, LINC00856, POLR2J3, TBCE. This implies these were problematic.
  # A more general approach for duplicated HGNCs:
  cat("  Handling duplicated HGNC symbols by keeping the one with highest average expression...\n")
  if (any(duplicated(FinalGeneFilter$hgnc_symbol))) {
    # Calculate row means for expression values (columns 3 to end)
    # Assuming sample columns start from column 3 after ensembl_gene_id and hgnc_symbol
    sample_cols_start_idx <- 3
    FinalGeneFilter$row_mean_expr <- rowMeans(FinalGeneFilter[, sample_cols_start_idx:ncol(FinalGeneFilter)], na.rm = TRUE)
    FinalGeneFilter <- FinalGeneFilter %>%
      dplyr::group_by(hgnc_symbol) %>%
      dplyr::filter(row_mean_expr == max(row_mean_expr)) %>%
      dplyr::ungroup() %>%
      dplyr::distinct(hgnc_symbol, .keep_all = TRUE) # Ensure unique HGNC
    FinalGeneFilter$row_mean_expr <- NULL # Remove temporary column
    cat("  Dimensions after de-duplicating HGNC symbols:", dim(FinalGeneFilter)[1], "genes.\n")
  }

  # RMD's specific removals (if the general deduplication isn't enough or different logic was intended)
  # FinalGeneFilter <- FinalGeneFilter[!(FinalGeneFilter$hgnc_symbol %in% c("STRA6LP", "LINC00856", "POLR2J3", "TBCE")),]

  # Final check: Ensure Ensembl IDs are unique as they should be primary keys from counts
  if (any(duplicated(FinalGeneFilter$ensembl_gene_id))) {
    cat("  Warning: Duplicated ensembl_gene_id found after filtering. Keeping first occurrence.\n")
    FinalGeneFilter <- FinalGeneFilter[!duplicated(FinalGeneFilter$ensembl_gene_id), ]
  }
  cat("  Dimensions after ensuring unique Ensembl IDs:", dim(FinalGeneFilter)[1], "genes.\n")

  # The RMD had a `keep = rowSums(FinalGeneFilter[2:13] >1) >= 6` line here.
  # This seems redundant if filtMddInMe was already filtered correctly.
  # The column indices also need adjustment if 'hgnc_symbol' is column 2.
  # If FinalGeneFilter has ensembl_gene_id, hgnc_symbol, then samples:
  # Sample columns start from index 3.
  # This filtering step might have been intended for after joining with gene symbols.
  cat("  Applying count filter again (rowSums > 1 in >= 6 samples) on the symbol-mapped data...\n")
  count_columns <- FinalGeneFilter[, 3:ncol(FinalGeneFilter)] # Assuming 1=ensembl_id, 2=hgnc_symbol
  keep_after_join <- rowSums(count_columns > 1) >= 6
  FinalGeneFilter <- FinalGeneFilter[keep_after_join, ]
  cat("  Final dimensions of 'FinalGeneFilter':", dim(FinalGeneFilter)[1], "genes,", dim(FinalGeneFilter)[2], "columns.\n")

  saveRDS(FinalGeneFilter, rds_GSE193417_filtered_symbols) # This is pre-column name fixing
  cat("  'FinalGeneFilter' (pre-column name fixing) saved to RDS.\n")

  # Fix sample column names (remove .bam suffix if present, etc.)
  # The RMD code:
  # testFinalGeneFilter <- FinalGeneFilter
  # for (i in 3:length(colnames(FinalGeneFilter))) {
  #   testnames <- strsplit((colnames(FinalGeneFilter)[i]), ".bam")
  #   colnames(testFinalGeneFilter)[i] <- testnames[[1]][1] # Access first element
  # }
  # A more robust way using gsub or stringr:
  testFinalGeneFilter <- FinalGeneFilter
  original_colnames <- colnames(testFinalGeneFilter)
  # Sample columns start from 3rd position
  new_colnames <- original_colnames
  for (i in 3:length(original_colnames)) {
    # Remove common suffixes or patterns. Example: ".bam", "_aligned.bam", "_counts.txt"
    # The RMD specifically targets ".bam". Let's make it more general or specific to GSE193417 names.
    # Sample names in GSE193417 raw counts file are like "CRH-Hu1001", "CRH-Hu103" etc.
    # It seems the RMD's `strsplit((colnames(FinalGeneFilter)[i]), ".bam")` was for a different upstream filename format.
    # For the current `filtMddInMe` from `readr::read_tsv`, colnames are already clean.
    # If they were like 'CRH-Hu1001.bam', this would be:
    # new_colnames[i] <- gsub("\\.bam$", "", original_colnames[i])
    # new_colnames[i] <- gsub("_.*", "", original_colnames[i]) # If more general cleaning needed
  }
  # If the column names from filtMddInMe are already clean, this step might not change anything.
  # Let's check the column names of FinalGeneFilter:
  cat("  Column names of FinalGeneFilter before potential renaming:\n")
  print(colnames(FinalGeneFilter))
  # The provided RMD's `strsplit` logic implies it expected sample names like `XXX.bam`.
  # The actual column names from the tsv are already clean (e.g., "CRH-Hu1001").
  # So, the renaming loop from RMD might not be necessary or might need adjustment
  # if the input `mddInMe` had different column names.
  # For now, assuming colnames are as from `read_tsv` of the NCBI file.
  # If they were `CRH-Hu1001.bamAligned.sortedByCoord.out.bam`, the RMD's split was too simple.
  # Let's assume the current `testFinalGeneFilter` colnames are what we need.

  colnames(testFinalGeneFilter) <- new_colnames
  cat("  Final column names for analysis ('testFinalGeneFilter'):\n")
  print(colnames(testFinalGeneFilter))

  saveRDS(testFinalGeneFilter, rds_GSE193417_finalfiltered_symbols)
  cat("  'testFinalGeneFilter' (with potentially fixed column names) saved to RDS.\n")
} else {
  cat("  Loading 'FinalGeneFilter' (pre-column name fixing) from RDS...\n")
  FinalGeneFilter <- readRDS(rds_GSE193417_filtered_symbols)
  cat("  Loading 'testFinalGeneFilter' (with potentially fixed column names) from RDS...\n")
  testFinalGeneFilter <- readRDS(rds_GSE193417_finalfiltered_symbols)
  cat("  Loaded 'testFinalGeneFilter' dimensions:", dim(testFinalGeneFilter)[1], "genes,", dim(testFinalGeneFilter)[2], "columns.\n")
}
cat("  Head of 'testFinalGeneFilter' (first 6 rows, first 5 columns):\n")
print(testFinalGeneFilter[1:6, 1:5])


# 1.6: Normalize counts using edgeR (TMM normalization)
cat("\nStep 1.6: Normalizing gene counts using edgeR (TMM method)...\n")
if (!file.exists(rds_GSE193417_normalized_datastruct)) {
  cat("  Creating DGEList object...\n")
  # Counts data are in columns 3 to end of testFinalGeneFilter
  # Ensure row names are set for DGEList if needed, e.g., by hgnc_symbol
  # The RMD used testFinalGeneFilter[3:14] which implies 12 sample columns.
  # Let's verify: ncol(testFinalGeneFilter) should be 2 (IDs) + 12 (samples) = 14
  num_sample_cols <- ncol(testFinalGeneFilter) - 2
  cat("  Number of sample columns detected:", num_sample_cols, "\n")

  # Create a DGEList object. Use ensembl_gene_id as row names for uniqueness,
  # but keep hgnc_symbol for later annotation.
  counts_for_dgelist <- testFinalGeneFilter[, 3:ncol(testFinalGeneFilter)]
  rownames(counts_for_dgelist) <- testFinalGeneFilter$ensembl_gene_id

  dge <- edgeR::DGEList(
    counts = counts_for_dgelist,
    genes = testFinalGeneFilter[, c("ensembl_gene_id", "hgnc_symbol")]
  )
  cat("  DGEList object created. Genes:", nrow(dge), "Samples:", ncol(dge), "\n")

  cat("  Calculating normalization factors (TMM)...\n")
  dge <- edgeR::calcNormFactors(dge, method = "TMM")
  cat("  Normalization factors:\n")
  print(dge$samples$norm.factors)

  cat("  Calculating Counts Per Million (CPM)...\n")
  normalized_counts_pmil <- edgeR::cpm(dge, log = FALSE) # Raw CPMs
  # If logCPM needed for MDS/density: log2(edgeR::cpm(dge, prior.count=2))
  cat("  Dimensions of normalized CPM matrix:", dim(normalized_counts_pmil)[1], "genes,", dim(normalized_counts_pmil)[2], "samples.\n")
  print(head(normalized_counts_pmil[, 1:4]))
  saveRDS(normalized_counts_pmil, rds_GSE193417_normalized_CPM)
  cat("  Normalized CPM matrix saved to RDS.\n")

  # Create the final normalized data structure (similar to RMD's normalized_datastruct)
  # This will have ensembl_gene_id, hgnc_symbol, then normalized counts
  # Ensure rownames of normalized_counts_pmil match dge$genes$ensembl_gene_id
  normalized_datastruct <- cbind(dge$genes, normalized_counts_pmil)
  cat("  Final normalized data structure ('normalized_datastruct') dimensions:", dim(normalized_datastruct)[1], "genes,", dim(normalized_datastruct)[2], "columns.\n")
  print(head(normalized_datastruct[, 1:5]))
  saveRDS(normalized_datastruct, rds_GSE193417_normalized_datastruct)
  cat("  Normalized data structure saved to RDS.\n")
} else {
  cat("  Loading normalized CPM matrix from RDS...\n")
  normalized_counts_pmil <- readRDS(rds_GSE193417_normalized_CPM)
  cat("  Loading normalized data structure from RDS...\n")
  normalized_datastruct <- readRDS(rds_GSE193417_normalized_datastruct)
  cat("  Loaded 'normalized_datastruct' dimensions:", dim(normalized_datastruct)[1], "genes,", dim(normalized_datastruct)[2], "columns.\n")

  # Recreate DGEList if needed for downstream (e.g. if voom is used later, though RMD uses limma on CPMs)
  counts_for_dgelist <- testFinalGeneFilter[, 3:ncol(testFinalGeneFilter)]
  rownames(counts_for_dgelist) <- testFinalGeneFilter$ensembl_gene_id
  dge <- edgeR::DGEList(
    counts = counts_for_dgelist,
    genes = testFinalGeneFilter[, c("ensembl_gene_id", "hgnc_symbol")]
  )
  dge <- edgeR::calcNormFactors(dge, method = "TMM") # Recalculate for consistency if needed
}

cat("--- COMPLETED SECTION 1 ---\n\n")

# ---
# SECTION 1B: ERROR CHECK AND SANITY PRINTOUTS
# ---
cat("--- STARTING SECTION 1B: ERROR CHECK AND SANITY PRINTOUTS ---\n")
if (!file.exists(rds_GSE193417_normalized_datastruct) || !exists("normalized_datastruct")) {
  stop("CRITICAL ERROR: Normalized data structure (normalized_datastruct.rds or object) not found. Please re-run Section 1.")
} else {
  cat("Normalized data structure ('normalized_datastruct') successfully loaded or created.\n")
  cat("  Number of genes:", nrow(normalized_datastruct), "\n")
  cat("  Number of columns:", ncol(normalized_datastruct), "\n")
  cat("  Column names:", colnames(normalized_datastruct), "\n")
  cat("  Summary of a few numeric columns (e.g., first sample column):\n")
  # Assuming 3rd column is the first sample
  if (ncol(normalized_datastruct) >= 3) print(summary(normalized_datastruct[[3]])) else print("Not enough columns for summary")
}
cat("--- COMPLETED SECTION 1B ---\n\n")


# ---
# SECTION 2: DIFFERENTIAL GENE EXPRESSION ANALYSIS
# ---
# This section constructs the design matrix, performs DGE analysis using limma,
# and visualizes results (MDS, Density, Volcano, MA, Heatmap).
# ---
cat("--- STARTING SECTION 2: DIFFERENTIAL GENE EXPRESSION ANALYSIS ---\n")

# 2.1: Construct the Design Matrix
cat("Step 2.1: Constructing the design matrix...\n")
# The RMD code for 'samples' and 'design2' construction:
# Samples names are CRH-Hu1001 sgACC_MDD, CRH-Hu103 sgACC_control, etc.
# The sample names in `normalized_datastruct` (derived from `testFinalGeneFilter`) should be used.
# Let's get sample names from `normalized_datastruct` (cols 3 onwards)
sample_column_names_from_norm_data <- colnames(normalized_datastruct)[3:ncol(normalized_datastruct)]

# Define groups based on sample names (MDD vs Control)
# This needs to be robust. The paper mentions 6 MDD and 6 Control.
# Assuming sample names contain "MDD" or "control" hints, or rely on pData from ExpressionSet.
# The RMD's `samples` object was constructed manually/semi-manually.
# Let's try to derive it from pData(GSE193417s) or infer from names if consistent.
# pData(GSE193417s)$title or pData(GSE193417s)$source_name_ch1 usually has group info.
# Example: pData(GSE193417s)$`disease state:ch1` might be "Major Depressive Disorder" or "control"
# Let's ensure the order of samples in pData matches our `normalized_datastruct` columns.
# The columns in `normalized_datastruct` (from `testFinalGeneFilter`) are:
# "CRH-Hu1001" "CRH-Hu103"  "CRH-Hu1047" "CRH-Hu1086" "CRH-Hu513" "CRH-Hu600"
# "CRH-Hu615"  "CRH-Hu789"  "CRH-Hu809"  "CRH-Hu852"  "CRH-Hu863" "CRH-Hu943"
# These names should match rownames of pData(GSE193417s) or part of geo_accession.
# The pData has 'title' like "CRH-Hu1001 sgACC_MDD".

# Let's use the RMD's manual definition of 'samples' as it's explicit about group mapping.
# Column names of `normalized_counts_pmil` are the sample identifiers.
sample_ids_for_design <- colnames(normalized_counts_pmil)

# Create the 'samples' metadata data frame as in RMD (manually curated or inferred)
# This is crucial for correct design matrix.
# RMD mapping:
# CRH-Hu1001 sgACC_MDD -> MDD
# CRH-Hu103 sgACC_control -> control
# ...
# Based on the RMD's `samples[1,]` and `samples[2,]`
sample_expression_type <- character(length(sample_ids_for_design))
sample_group <- character(length(sample_ids_for_design))

# This mapping needs to be carefully established based on paper or metadata
# For GSE193417, the `title` in pData(GSE193417s) indicates group.
# E.g., "CRH-Hu1001 sgACC_MDD"
pheno_data <- pData(GSE193417s)
# Match sample_ids_for_design with pheno_data rownames or titles to get groups
# The column names of normalized_counts_pmil are `CRH-HuXXXX`
# The rownames of pheno_data are GEO sample IDs (GSMXXXXXX)
# The `title` column in pheno_data seems to be `CRH-HuXXXX sgACC_GROUP`

# Recreate 'design2' from the RMD (sample_expression, sample_group)
# The RMD's 'design2' had columns 'sample_expression' and 'sample_group'.
# 'sample_expression' was 'sgACC_MDD' or 'sgACC_control'.
# 'sample_group' was 'MDD' or 'control'.

# We need to map our `sample_ids_for_design` to these groups.
# Let's use the explicit mapping from the RMD's `samples` object logic.
design_df_data <- list()
for (id in sample_ids_for_design) {
  title_match <- pheno_data$title[grepl(id, pheno_data$title)]
  if (length(title_match) == 1) {
    if (grepl("MDD", title_match)) {
      design_df_data[[id]] <- data.frame(
        sample_id = id,
        sample_expression = "sgACC_MDD",
        sample_group = "MDD",
        stringsAsFactors = FALSE
      )
    } else if (grepl("control", title_match, ignore.case = TRUE)) {
      design_df_data[[id]] <- data.frame(
        sample_id = id,
        sample_expression = "sgACC_control",
        sample_group = "control",
        stringsAsFactors = FALSE
      )
    } else {
      warning(paste("Could not determine group for sample:", id, "from title:", title_match))
      design_df_data[[id]] <- data.frame(sample_id = id, sample_expression = NA, sample_group = NA)
    }
  } else {
    warning(paste("Could not uniquely map sample ID:", id, "to phenoData titles."))
    # Fallback to RMD's direct assignment structure if titles are ambiguous
    # This requires knowing the exact order and assignments.
    # For robustness, let's use the pattern from RMD's `samples[1,]` and `samples[2,]`
    # CRH-Hu1001 sgACC_MDD, CRH-Hu103 sgACC_control, CRH-Hu1047 sgACC_control, CRH-Hu1086 sgACC_control,
    # CRH-Hu513 sgACC_MDD, CRH-Hu600 sgACC_MDD, CRH-Hu615 sgACC_control, CRH-Hu789 sgACC_control,
    # CRH-Hu809 sgACC_MDD, CRH-Hu852 sgACC_control, CRH-Hu863 sgACC_MDD, CRH-Hu943 sgACC_MDD
    # This order corresponds to columns in `normalized_counts_pmil`

    # RMD's `samples[1,]` and `samples[2,]`
    rmd_sample_expressions <- c("sgACC_MDD", "sgACC_control", "sgACC_control", "sgACC_control", "sgACC_MDD", "sgACC_MDD", "sgACC_control", "sgACC_control", "sgACC_MDD", "sgACC_control", "sgACC_MDD", "sgACC_MDD")
    rmd_sample_groups <- c("MDD", "control", "control", "control", "MDD", "MDD", "control", "control", "MDD", "control", "MDD", "MDD")

    idx_in_rmd_order <- match(id, sample_ids_for_design) # This assumes sample_ids_for_design is in the RMD's implicit order
    if (!is.na(idx_in_rmd_order)) {
      design_df_data[[id]] <- data.frame(
        sample_id = id,
        sample_expression = rmd_sample_expressions[idx_in_rmd_order],
        sample_group = rmd_sample_groups[idx_in_rmd_order],
        stringsAsFactors = FALSE
      )
    } else {
      # This should not happen if sample_ids_for_design is correct
      design_df_data[[id]] <- data.frame(sample_id = id, sample_expression = NA, sample_group = NA)
    }
  }
}
design_df <- do.call(rbind, design_df_data)
rownames(design_df) <- design_df$sample_id

# Ensure the order of design_df matches columns of normalized_counts_pmil
design_df <- design_df[match(colnames(normalized_counts_pmil), design_df$sample_id), ]

cat("  Constructed design_df (equivalent to RMD's design2):\n")
print(design_df)

# Create the model matrix for limma. We are interested in the difference between MDD and control.
# So, sample_group is the factor of interest.
# It must be a factor for model.matrix.
design_df$sample_group <- factor(design_df$sample_group, levels = c("control", "MDD")) # Set control as reference

# The RMD used `~ design2$sample_expression`.
# If `sample_expression` is `sgACC_MDD` vs `sgACC_control`, this is also valid.
# Let's stick to `sample_group` as it's simpler for pairwise comparison.
# If `~design_df$sample_expression` was used, ensure its levels are set.
design_df$sample_expression <- factor(design_df$sample_expression, levels = c("sgACC_control", "sgACC_MDD"))

# Using sample_group as it's clearer for 'MDD vs control'
# design_model <- model.matrix(~ sample_group, data = design_df)
# Using sample_expression as in the RMD
design_model <- model.matrix(~sample_expression, data = design_df)

rownames(design_model) <- design_df$sample_id
cat("  Final design_model for limma:\n")
print(design_model)

# 2.2: Data visualization (MDS Plot, Density Curve)
cat("\nStep 2.2: Data visualization (MDS, Density)...\n")

# MDS Plot on log2 CPM values
# The RMD used `d <- testFinalGeneFilter[3:14]`. This is raw counts.
# For MDS, it's common to use log-transformed normalized counts.
# Let's use log2(normalized_counts_pmil + prior.count) for MDS
# A small prior.count (e.g., 0.25 or 1) avoids log(0)
log2_cpm_for_mds <- log2(normalized_counts_pmil + 1)

cat("  Generating MDS plot...\n")
# The groups for coloring: use 'sample_group' from design_df
mds_groups <- design_df$sample_group
# plotMDS function from limma or edgeR. edgeR::plotMDS is common for DGEList.
# limma::plotMDS works on a matrix.
# Ensure colors are consistent.
group_colors <- brewer.pal(length(levels(mds_groups)), "Set1")
names(group_colors) <- levels(mds_groups)

# Create a new plot window
if (.Platform$OS.type != "windows") X11() else windows()
limma::plotMDS(log2_cpm_for_mds,
  labels = colnames(log2_cpm_for_mds),
  col = group_colors[mds_groups],
  main = "MDS Plot of Normalized log2(CPM+1) Values for GSE193417"
)
legend("topleft", legend = levels(mds_groups), fill = group_colors, cex = 0.8)
cat("  MDS plot generated. Check graphics device.\n")


# Density Curve for Normalized Distribution Values (log2 CPM)
cat("  Generating density plot of log2(CPM+1) values...\n")
# Data for plot: log2_cpm_for_mds
ctsDenct <- apply(log2_cpm_for_mds, 2, density)
sample_names_for_plot <- colnames(log2_cpm_for_mds)

xlim_density <- range(sapply(ctsDenct, function(d) range(d$x)))
ylim_density <- range(sapply(ctsDenct, function(d) range(d$y)))

density_plot_colors <- rainbow(length(ctsDenct)) # One color per sample
density_plot_ltys <- rep(1, length(ctsDenct))

if (.Platform$OS.type != "windows") X11() else windows()
plot(0, 0,
  type = "n", # Empty plot
  xlim = xlim_density, ylim = ylim_density,
  xlab = "log2(CPM+1)", ylab = "Density",
  main = "Density Curves of Normalized log2(CPM+1) Values for GSE193417",
  cex.lab = 0.85
)

for (i in 1:length(ctsDenct)) {
  lines(ctsDenct[[i]], col = density_plot_colors[i], lty = density_plot_ltys[i])
}
legend("topright",
  legend = sample_names_for_plot, col = density_plot_colors,
  lty = density_plot_ltys, cex = 0.7, ncol = ifelse(length(ctsDenct) > 6, 2, 1),
  title = "Samples"
)
cat("  Density plot generated. Check graphics device.\n")


# 2.3: Perform Differential Gene Expression with limma
cat("\nStep 2.3: Performing DGE analysis with limma...\n")
# limma typically works well with log-transformed data (like logCPM from voom, or microarray data).
# The RMD used `ExpressionSet(normalized_counts_pmil)`. These are CPMs, not logCPMs.
# While limma can run on CPMs, it's more standard to use logCPMs for linear modeling.
# If `voom` was used on raw counts (in `dge` object), it would produce logCPMs with precision weights.
# Let's use log2(normalized_counts_pmil + prior.count) for the linear model fit, consistent with MDS/Density.
# prior.count value for log transformation can be small, e.g., 0.25 or 1. Let's use 1 for simplicity.
log2_cpm_for_limma <- log2(normalized_counts_pmil + 1)

cat("  Creating ExpressionSet object with log2(CPM+1)...\n")
# Ensure phenoData for ExpressionSet matches sample order
# Rownames of phenoData should be sample names (colnames of expression matrix)
# design_df already has sample_id as rownames and is ordered.
phenoData_for_eset <- AnnotatedDataFrame(data = design_df)
minimal_set <- ExpressionSet(
  assayData = log2_cpm_for_limma,
  phenoData = phenoData_for_eset
)
cat("  ExpressionSet created.\n")
print(minimal_set)

cat("  Fitting linear model...\n")
fit <- limma::lmFit(minimal_set, design_model)
cat("  Linear model fit completed. Coefficients found:\n")
print(colnames(fit$coefficients))

cat("  Applying empirical Bayes smoothing...\n")
# The RMD used trend=TRUE. This is generally good if there's a mean-variance trend.
fit2 <- limma::eBayes(fit, trend = TRUE, robust = TRUE) # Added robust=TRUE
cat("  Empirical Bayes smoothing completed.\n")

# Extract top differentially expressed genes
# The coefficient of interest depends on how `design_model` was set up.
# design_model <- model.matrix(~ sample_expression, data = design_df)
# sample_expression levels: "sgACC_control", "sgACC_MDD"
# The second coefficient `sample_expression_sgACC_MDD` represents MDD vs Control.
coef_of_interest <- colnames(design_model)[ncol(design_model)] # Usually the last one for simple designs
cat("  Coefficient of interest for DGE:", coef_of_interest, "\n")

# Get all genes, adjust p-values using Benjamini-Hochberg
topfit <- limma::topTable(fit2,
  coef = coef_of_interest,
  adjust.method = "BH",
  number = Inf
) # Inf gets all genes
cat("  Extracted DGE results for", nrow(topfit), "genes.\n")
cat("  Head of DGE results ('topfit'):\n")
print(head(topfit))

# Add gene symbols to the results
# `topfit` rownames are Ensembl IDs (from rownames of log2_cpm_for_limma)
# `dge$genes` has ensembl_gene_id and hgnc_symbol
# Create a mapping data frame from dge$genes
gene_annotation_df <- dge$genes[, c("ensembl_gene_id", "hgnc_symbol")]
# Merge with topfit. topfit rownames are ensembl_gene_id.
output_hits <- merge(gene_annotation_df, topfit, by.x = "ensembl_gene_id", by.y = "row.names", all.y = TRUE)
# Rename first column if it becomes "Row.names" to something like "ID" or "ensembl_gene_id"
# The by.y="row.names" makes rownames a column, so by.x="ensembl_gene_id" adds it.
# Let's check output_hits structure
# It will have 'ensembl_gene_id' from gene_annotation_df, and then other columns from topfit.
# If gene_annotation_df$ensembl_gene_id was not unique, this could be an issue. It should be.
# Ensure order from topfit is preserved if needed, or re-sort.
# RMD sorts by unadjusted P.Value later.

# Rename columns for clarity if needed e.g. P.Value, adj.P.Val, logFC
# output_hits already has standard limma names: logFC, AveExpr, t, P.Value, adj.P.Val, B
cat("  Merged DGE results with gene symbols ('output_hits').\n")
cat("  Head of 'output_hits':\n")
print(head(output_hits))
cat("  Number of NAs in P.Value:", sum(is.na(output_hits$P.Value)), "\n")
cat("  Number of NAs in adj.P.Val:", sum(is.na(output_hits$adj.P.Val)), "\n")

# Sort by unadjusted P.Value (as in RMD)
output_hits <- output_hits[order(output_hits$P.Value, decreasing = FALSE), ]
cat("  'output_hits' sorted by P.Value. Head:\n")
print(head(output_hits))
cat("  Number of genes with P.Value < 0.05:", sum(output_hits$P.Value < 0.05, na.rm = TRUE), "\n")
cat("  Number of genes with adj.P.Val < 0.05 (FDR < 5%):", sum(output_hits$adj.P.Val < 0.05, na.rm = TRUE), "\n")
# The RMD's question used P.Value < 0.05 for "significantly DE" before correction.
# For BH correction, it's `adj.P.Val`.


# 2.4: Visualize DGE Results
cat("\nStep 2.4: Visualizing DGE results (Volcano, MA, Heatmap)...\n")

# Volcano Plot
cat("  Generating Volcano Plot...\n")
# Add a column for differential expression status
vol_plot_data <- output_hits
vol_plot_data$diffexpressed <- "Not Significant"
# Thresholds: logFC > 1 (or < -1) AND P.Value < 0.05 (RMD used P.Value, often adj.P.Val is used)
# Let's use adj.P.Val for significance, and a logFC threshold.
LFC_threshold <- 1
P_adj_threshold <- 0.05 # For stricter significance

vol_plot_data$diffexpressed[which(vol_plot_data$logFC > LFC_threshold & vol_plot_data$adj.P.Val < P_adj_threshold)] <- "Up Regulated"
vol_plot_data$diffexpressed[which(vol_plot_data$logFC < -LFC_threshold & vol_plot_data$adj.P.Val < P_adj_threshold)] <- "Down Regulated"
cat("  Differential expression status distribution:\n")
print(table(vol_plot_data$diffexpressed))

# Add labels for top genes
# Create a 'delabel' column as in RMD
vol_plot_data$delabel <- NA
# Label genes that are significant (Up or Down Regulated)
genes_to_label_idx <- which(vol_plot_data$diffexpressed != "Not Significant")
vol_plot_data$delabel[genes_to_label_idx] <- vol_plot_data$hgnc_symbol[genes_to_label_idx]
# Optionally, label only the very top N genes to avoid clutter
# For example, top 10 by smallest adj.P.Val among significant ones
# top_n_to_label <- 10
# significant_genes_for_labeling <- vol_plot_data[genes_to_label_idx,]
# significant_genes_for_labeling <- significant_genes_for_labeling[order(significant_genes_for_labeling$adj.P.Val),]
# vol_plot_data$delabel <- NA # Reset
# if(nrow(significant_genes_for_labeling) > 0) {
#   label_indices_in_vol_plot_data <- match(head(significant_genes_for_labeling$ensembl_gene_id, top_n_to_label), vol_plot_data$ensembl_gene_id)
#   vol_plot_data$delabel[label_indices_in_vol_plot_data] <- vol_plot_data$hgnc_symbol[label_indices_in_vol_plot_data]
# }


# RMD used P.Value < 0.05 for its example genes (MT-CO1, GAS6 etc.)
# Let's make a version with P.Value for coloring like RMD, and another with adj.P.Val
# RMD's 'vol_plot' definition:
# vol_plot$diffexpressed[vol_plot$logFC > 1 & vol_plot$P.Value < 0.05] <- "Up Regulated"
# vol_plot$diffexpressed[vol_plot$logFC < -1 & vol_plot$P.Value < 0.05] <- "Down Regulated"
# (Note: RMD had < 1 for downregulated, should be < -1 for magnitude)
# This is for consistency with RMD's example text for highlighted genes
vol_plot_data_rmd_criteria <- output_hits
vol_plot_data_rmd_criteria$diffexpressed <- "Not Significant"
vol_plot_data_rmd_criteria$diffexpressed[vol_plot_data_rmd_criteria$logFC > 1 & vol_plot_data_rmd_criteria$P.Value < 0.05] <- "Up Regulated"
vol_plot_data_rmd_criteria$diffexpressed[vol_plot_data_rmd_criteria$logFC < -1 & vol_plot_data_rmd_criteria$P.Value < 0.05] <- "Down Regulated" # Corrected from RMD's <1
vol_plot_data_rmd_criteria$delabel <- NA
vol_plot_data_rmd_criteria$delabel[vol_plot_data_rmd_criteria$diffexpressed != "Not Significant"] <- vol_plot_data_rmd_criteria$hgnc_symbol[vol_plot_data_rmd_criteria$diffexpressed != "Not Significant"]


if (.Platform$OS.type != "windows") X11() else windows()
p_volcano <- ggplot(
  data = vol_plot_data_rmd_criteria,
  aes(
    x = logFC, y = -log10(P.Value),
    col = diffexpressed, label = delabel
  )
) +
  geom_point(alpha = 0.6, size = 1.5) +
  theme_minimal(base_size = 12) +
  ggrepel::geom_text_repel(
    max.overlaps = 15, size = 3,
    box.padding = 0.5, point.padding = 0.2,
    min.segment.length = 0
  ) + # Added min.segment.length
  scale_color_manual(
    values = c("blue", "black", "red"), # Down, Not Sig, Up
    name = "Differential Expression",
    labels = c("Down Regulated", "Not Significant", "Up Regulated")
  ) +
  geom_vline(xintercept = c(-LFC_threshold, LFC_threshold), col = "grey50", linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), col = "grey50", linetype = "dashed") +
  labs(
    title = "Volcano Plot: Differentially Expressed Genes (GSE193417)",
    subtitle = paste0("logFC threshold: +/-", LFC_threshold, ", P.Value threshold: 0.05 (unadjusted)"),
    x = "log2 Fold Change (MDD vs Control)",
    y = "-log10(P.Value)"
  ) +
  theme(legend.position = "bottom")
print(p_volcano)
cat("  Volcano plot generated. Check graphics device.\n")


# MA Plot
cat("  Generating MA Plot...\n")
# limma::plotMA uses the 'fit2' object.
# Status can highlight genes based on a condition.
# RMD: status=output_hits$P.Value < 0.05
# Let's use the same criteria for consistency.
# Ensure output_hits is ordered same as genes in fit2 (by ensembl_gene_id)
# The 'fit2' object has genes in the order of the input expression matrix.
# `output_hits` was sorted by P.Value. We need to match its P.Values back to `fit2$genes` order.
# `fit2$genes` is not populated by default unless `lmFit` was given `dge$genes`.
# `minimal_set` was created from `log2_cpm_for_limma` which had ensembl_gene_ids as rownames.
# `fit2` will have these as `fit2$Amean` names or similar.
# The `plotMA` function can take an MArrayLM object directly.
status_for_maplot <- rep(FALSE, nrow(fit2))
# Get ensembl_ids from fit2 (usually from rownames of data that went into lmFit)
ensembl_ids_in_fit2 <- rownames(log2_cpm_for_limma) # Or featureNames(minimal_set)

# Match these with output_hits to get P.Values
p_values_ordered_as_fit2 <- output_hits$P.Value[match(ensembl_ids_in_fit2, output_hits$ensembl_gene_id)]
status_for_maplot_condition <- p_values_ordered_as_fit2 < 0.05 & !is.na(p_values_ordered_as_fit2)

if (.Platform$OS.type != "windows") X11() else windows()
limma::plotMA(fit2,
  coef = coef_of_interest,
  status = status_for_maplot_condition,
  values = c("TRUE"), # Value in status that indicates significance
  col = c("red"), pch = c(16), # Color and point type for significant
  main = "MA Plot: Differentially Expressed Genes (GSE193417)",
  xlab = "Average log2 Expression", ylab = "log2 Fold Change (MDD vs Control)",
  legend = "bottomright"
) # Added legend position
abline(h = c(-LFC_threshold, LFC_threshold), col = "dodgerblue", lty = 2)
cat("  MA plot generated. Check graphics device.\n")


# Heatmap of top differentially expressed genes
cat("  Generating Heatmap for top DE genes...\n")
# RMD: topHitHM <- output_hits$Row.names[output_hits$P.Value < 0.05]
# `Row.names` is not a column in `output_hits` after merge if done correctly. It uses `ensembl_gene_id`.
# Let's select top N genes by adj.P.Val or P.Value.
# Using P.Value < 0.05 as in RMD.
# `output_hits` is already sorted by P.Value.
top_de_genes_for_heatmap_ensids <- output_hits$ensembl_gene_id[output_hits$P.Value < 0.05 & !is.na(output_hits$P.Value)]
cat("  Number of genes for heatmap (P.Value < 0.05):", length(top_de_genes_for_heatmap_ensids), "\n")

# Select a subset if too many, e.g., top 50
if (length(top_de_genes_for_heatmap_ensids) > 50) {
  cat("  More than 50 genes, selecting top 50 by P.Value for heatmap.\n")
  top_de_genes_for_heatmap_ensids <- head(top_de_genes_for_heatmap_ensids, 50)
} else if (length(top_de_genes_for_heatmap_ensids) == 0) {
  cat("  No genes found with P.Value < 0.05. Skipping heatmap.\n")
}

if (length(top_de_genes_for_heatmap_ensids) > 0) {
  # Get normalized expression data for these genes (log2 CPM)
  # `normalized_counts_pmil` has ensembl_gene_ids as rownames.
  heatmap_matrix_data_raw <- log2_cpm_for_limma[rownames(log2_cpm_for_limma) %in% top_de_genes_for_heatmap_ensids, ]

  # Add HGNC symbols as rownames for heatmap display
  hgnc_for_heatmap_genes <- output_hits$hgnc_symbol[match(rownames(heatmap_matrix_data_raw), output_hits$ensembl_gene_id)]
  rownames(heatmap_matrix_data_raw) <- make.unique(hgnc_for_heatmap_genes) # Ensure unique rownames

  # Scale data by row (z-score for each gene across samples)
  heatmap_matrix_scaled <- t(scale(t(heatmap_matrix_data_raw)))
  # Remove rows with NaN if scaling produced them (e.g. zero variance genes, though unlikely for DEGs)
  heatmap_matrix_scaled <- heatmap_matrix_scaled[complete.cases(heatmap_matrix_scaled), ]

  cat("  Dimensions of matrix for heatmap:", dim(heatmap_matrix_scaled)[1], "genes,", dim(heatmap_matrix_scaled)[2], "samples.\n")

  # Define color ramp for heatmap
  heatmap_colors <- circlize::colorRamp2(
    breaks = c(min(heatmap_matrix_scaled, na.rm = T), 0, max(heatmap_matrix_scaled, na.rm = T)),
    colors = c("blue", "white", "red")
  )

  # Create column annotation for heatmap (sample groups)
  col_annotation_df <- data.frame(Group = design_df$sample_group)
  rownames(col_annotation_df) <- colnames(heatmap_matrix_scaled)
  col_annotation_colors <- list(Group = c("control" = "skyblue", "MDD" = "tomato"))

  heatmap_col_ann <- ComplexHeatmap::HeatmapAnnotation(
    df = col_annotation_df,
    col = col_annotation_colors,
    annotation_name_side = "left"
  )

  if (.Platform$OS.type != "windows") X11() else windows()
  complex_heatmap <- ComplexHeatmap::Heatmap(
    heatmap_matrix_scaled,
    name = "Z-score\n(log2 CPM)",
    col = heatmap_colors,
    cluster_rows = TRUE, show_row_dend = TRUE, row_dend_reorder = TRUE,
    cluster_columns = TRUE, show_column_dend = TRUE, column_dend_reorder = TRUE,
    show_row_names = ifelse(nrow(heatmap_matrix_scaled) <= 50, TRUE, FALSE), # Show row names if not too many
    row_names_gp = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 10),
    top_annotation = heatmap_col_ann,
    column_title = "Heatmap of Top Differentially Expressed Genes (P.Value < 0.05)",
    column_title_gp = gpar(fontsize = 12, fontface = "bold"),
    show_heatmap_legend = TRUE
  )
  draw(complex_heatmap, heatmap_legend_side = "right", annotation_legend_side = "right")
  cat("  Heatmap generated. Check graphics device.\n")
}
cat("--- COMPLETED SECTION 2 ---\n\n")


# ---
# SECTION 3: ANSWERS TO DGE QUESTIONS (Programmatic parts)
# ---
cat("--- STARTING SECTION 3: ANSWERS TO DGE QUESTIONS ---\n")

# Q1a: P-values for each gene are in `output_hits$P.Value`.
cat("Q1a: P-values for each gene are available in the 'output_hits' data frame, column 'P.Value'.\n")
cat("     Example (first 6 genes from sorted list):\n")
print(head(data.frame(
  ensembl_id = output_hits$ensembl_gene_id,
  hgnc = output_hits$hgnc_symbol,
  P.Value = output_hits$P.Value
)))

# Q1b: How many genes were significantly differentially expressed? (Using P.Value < 0.05 as per RMD)
num_sig_de_p05 <- sum(output_hits$P.Value < 0.05, na.rm = TRUE)
total_genes_tested <- nrow(output_hits)
percentage_sig_de_p05 <- (num_sig_de_p05 / total_genes_tested) * 100
cat(paste0("\nQ1b: Number of genes significantly DE (P.Value < 0.05, unadjusted): ", num_sig_de_p05, "\n"))
cat(paste0("     This is ", round(percentage_sig_de_p05, 2), "% of the ", total_genes_tested, " tested genes.\n"))

# Q1c: Thresholds used: P.Value < 0.05 (unadjusted).
cat("\nQ1c: Threshold used for 'significantly DE' in Q1b was P.Value < 0.05 (unadjusted), consistent with RMD's initial exploration.\n")
cat("     For more stringent analysis, adj.P.Val (FDR) is typically used (see Q2b).\n")


# Q2a: Multiple hypothesis correction method: Benjamini-Hochberg (BH) was used in topTable.
cat("\nQ2a: Multiple hypothesis correction method used was Benjamini-Hochberg (BH).\n")
cat("     This was specified via adjust.method='BH' in limma::topTable.\n")
cat("     BH method controls the False Discovery Rate (FDR).\n")

# Q2b: How many genes passed correction? (Typically adj.P.Val < threshold, e.g., 0.05 or 0.10)
fdr_threshold <- 0.05
num_genes_passed_fdr <- sum(output_hits$adj.P.Val < fdr_threshold, na.rm = TRUE)
cat(paste0("\nQ2b: Number of genes passing BH correction (adj.P.Val < ", fdr_threshold, "): ", num_genes_passed_fdr, "\n"))
# The RMD stated 0 genes passed correction. This depends on the data and chosen FDR. Let's see.

cat("--- COMPLETED SECTION 3 ---\n\n")


# ---
# SECTION 4: THRESHOLDED OVER-REPRESENTATION ANALYSIS (ORA)
# ---
# This section prepares gene lists for ORA and uses g:Profiler.
# ---
cat("--- STARTING SECTION 4: THRESHOLDED ORA ---\n")

# 4.1: Create lists of all, upregulated, and downregulated genes for ORA
cat("Step 4.1: Creating gene lists for ORA...\n")
# Using adj.P.Val < 0.05 and logFC threshold for defining DE genes for ORA
# This is a common practice for ORA to reduce noise.
ORA_LFC_THRESHOLD <- LFC_threshold # Use same LFC as volcano for consistency, or set new one
ORA_ADJP_THRESHOLD <- 0.05

# All significantly DE genes (using HGNC symbols for g:Profiler query)
all_de_genes_ora_hgnc <- output_hits$hgnc_symbol[
  which(output_hits$adj.P.Val < ORA_ADJP_THRESHOLD & abs(output_hits$logFC) > ORA_LFC_THRESHOLD & !is.na(output_hits$hgnc_symbol))
]
all_de_genes_ora_hgnc <- unique(all_de_genes_ora_hgnc[!is.na(all_de_genes_ora_hgnc) & all_de_genes_ora_hgnc != ""])
cat(paste0("  Number of all significantly DE genes for ORA (adj.P.Val < ", ORA_ADJP_THRESHOLD, ", |logFC| > ", ORA_LFC_THRESHOLD, "): ", length(all_de_genes_ora_hgnc), "\n"))

# Upregulated genes
up_regulated_ora_hgnc <- output_hits$hgnc_symbol[
  which(output_hits$adj.P.Val < ORA_ADJP_THRESHOLD & output_hits$logFC > ORA_LFC_THRESHOLD & !is.na(output_hits$hgnc_symbol))
]
up_regulated_ora_hgnc <- unique(up_regulated_ora_hgnc[!is.na(up_regulated_ora_hgnc) & up_regulated_ora_hgnc != ""])
cat(paste0("  Number of upregulated DE genes for ORA: ", length(up_regulated_ora_hgnc), "\n"))

# Downregulated genes
down_regulated_ora_hgnc <- output_hits$hgnc_symbol[
  which(output_hits$adj.P.Val < ORA_ADJP_THRESHOLD & output_hits$logFC < -ORA_LFC_THRESHOLD & !is.na(output_hits$hgnc_symbol))
]
down_regulated_ora_hgnc <- unique(down_regulated_ora_hgnc[!is.na(down_regulated_ora_hgnc) & down_regulated_ora_hgnc != ""])
cat(paste0("  Number of downregulated DE genes for ORA: ", length(down_regulated_ora_hgnc), "\n"))

# Background gene set for ORA (all genes tested in DGE analysis with valid HGNC symbols)
background_genes_hgnc <- unique(output_hits$hgnc_symbol[!is.na(output_hits$hgnc_symbol) & output_hits$hgnc_symbol != ""])
cat(paste0("  Number of background genes for ORA (universe): ", length(background_genes_hgnc), "\n"))

# Save gene lists to files (as in RMD, though gprofiler2 can take vectors directly)
# Using Ensembl IDs as RMD did for text files.
all_up_ens_ora <- output_hits$ensembl_gene_id[which(output_hits$adj.P.Val < ORA_ADJP_THRESHOLD & output_hits$logFC > ORA_LFC_THRESHOLD)]
all_down_ens_ora <- output_hits$ensembl_gene_id[which(output_hits$adj.P.Val < ORA_ADJP_THRESHOLD & output_hits$logFC < -ORA_LFC_THRESHOLD)]

# Writing to files (optional, gprofiler2 R package takes vectors)
# write.table(all_up_ens_ora, file = file.path(data_dir, "ora_up_regulated_ens.txt"), sep="\n", row.names = FALSE, col.names = FALSE, quote = FALSE)
# write.table(up_regulated_ora_hgnc, file = file.path(data_dir, "ora_up_regulated_hgnc.txt"), sep="\n", row.names = FALSE, col.names = FALSE, quote = FALSE)
# write.table(all_down_ens_ora, file = file.path(data_dir, "ora_down_regulated_ens.txt"), sep="\n", row.names = FALSE, col.names = FALSE, quote = FALSE)
# write.table(down_regulated_ora_hgnc, file = file.path(data_dir, "ora_down_regulated_hgnc.txt"), sep="\n", row.names = FALSE, col.names = FALSE, quote = FALSE)


# 4.2: Perform ORA using g:Profiler2
cat("\nStep 4.2: Performing ORA with g:Profiler2...\n")
# Define common sources for g:Profiler
gprofiler_sources <- c("GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC", "HPA", "HP") # WP was in RMD, HP is Human Phenotype

# ORA for all DE genes
cat("  Running ORA for all significantly DE genes...\n")
if (length(all_de_genes_ora_hgnc) > 0) {
  gostres_all_de <- gost(
    query = all_de_genes_ora_hgnc,
    organism = "hsapiens",
    ordered_query = FALSE,
    multi_query = FALSE,
    significant = TRUE, # Keep only significant results
    exclude_iea = FALSE, # Include electronic annotations
    measure_underrepresentation = FALSE,
    evcodes = TRUE, # Return evidence codes
    user_threshold = 0.05, # P-value threshold for terms
    correction_method = "g_SCS", # g:Profiler's default
    domain_scope = "annotated",
    custom_bg = background_genes_hgnc,
    sources = gprofiler_sources
  )

  cat("  ORA for all DE genes completed.\n")
  if (!is.null(gostres_all_de) && nrow(gostres_all_de$result) > 0) {
    cat("  Significant terms found for all DE genes:\n")
    print(head(gostres_all_de$result[, c("term_id", "term_name", "p_value", "source")]))
    if (.Platform$OS.type != "windows") X11() else windows()
    p_gost_all <- gostplot(gostres_all_de, capped = TRUE, interactive = FALSE) # interactive=FALSE for script
    print(p_gost_all)
    title(main = "g:Profiler ORA Results: All DE Genes", outer = TRUE, line = -1)
    cat("  g:Profiler plot for all DE genes generated.\n")
  } else {
    cat("  No significant terms found for all DE genes.\n")
  }
} else {
  cat("  Skipping ORA for all DE genes: No genes in the list.\n")
}

# ORA for upregulated DE genes
cat("\n  Running ORA for upregulated DE genes...\n")
if (length(up_regulated_ora_hgnc) > 0) {
  gostres_up_de <- gost(
    query = up_regulated_ora_hgnc,
    organism = "hsapiens",
    custom_bg = background_genes_hgnc,
    sources = gprofiler_sources,
    user_threshold = 0.05, correction_method = "g_SCS"
  )

  cat("  ORA for upregulated DE genes completed.\n")
  if (!is.null(gostres_up_de) && nrow(gostres_up_de$result) > 0) {
    cat("  Significant terms found for upregulated DE genes:\n")
    print(head(gostres_up_de$result[, c("term_id", "term_name", "p_value", "source")]))
    if (.Platform$OS.type != "windows") X11() else windows()
    p_gost_up <- gostplot(gostres_up_de, capped = TRUE, interactive = FALSE)
    print(p_gost_up)
    title(main = "g:Profiler ORA Results: Upregulated DE Genes", outer = TRUE, line = -1)
    cat("  g:Profiler plot for upregulated DE genes generated.\n")
  } else {
    cat("  No significant terms found for upregulated DE genes.\n")
  }
} else {
  cat("  Skipping ORA for upregulated DE genes: No genes in the list.\n")
}

# ORA for downregulated DE genes
cat("\n  Running ORA for downregulated DE genes...\n")
if (length(down_regulated_ora_hgnc) > 0) {
  gostres_down_de <- gost(
    query = down_regulated_ora_hgnc,
    organism = "hsapiens",
    custom_bg = background_genes_hgnc,
    sources = gprofiler_sources,
    user_threshold = 0.05, correction_method = "g_SCS"
  )

  cat("  ORA for downregulated DE genes completed.\n")
  if (!is.null(gostres_down_de) && nrow(gostres_down_de$result) > 0) {
    cat("  Significant terms found for downregulated DE genes:\n")
    print(head(gostres_down_de$result[, c("term_id", "term_name", "p_value", "source")]))
    if (.Platform$OS.type != "windows") X11() else windows()
    p_gost_down <- gostplot(gostres_down_de, capped = TRUE, interactive = FALSE)
    print(p_gost_down)
    title(main = "g:Profiler ORA Results: Downregulated DE Genes", outer = TRUE, line = -1)
    cat("  g:Profiler plot for downregulated DE genes generated.\n")
  } else {
    cat("  No significant terms found for downregulated DE genes.\n")
  }
} else {
  cat("  Skipping ORA for downregulated DE genes: No genes in the list.\n")
}

cat("--- COMPLETED SECTION 4 ---\n\n")


# ---
# SECTION 5: ANSWERS TO ORA QUESTIONS (Programmatic parts)
# ---
cat("--- STARTING SECTION 5: ANSWERS TO ORA QUESTIONS ---\n")

# Q_ORA_1: Method chosen: Over-Representation Analysis (ORA) using g:Profiler2.
cat("Q_ORA_1: Method chosen: Over-Representation Analysis (ORA) via g:Profiler2 R package.\n")
cat("         Why: g:Profiler is comprehensive, supports many databases, handles multiple testing correction, and allows custom backgrounds.\n")

# Q_ORA_2a: Annotation data used: GO (BP, MF, CC), KEGG, REAC, HPA, HP.
cat("\nQ_ORA_2a: Annotation data sources used: GO (Biological Process, Molecular Function, Cellular Component), KEGG, Reactome, Human Protein Atlas, Human Phenotype Ontology.\n")
cat("           Why: These cover a broad range of biological functions, pathways, and cellular locations relevant to gene expression studies.\n")

# Q_ORA_2b: Version of annotation: g:Profiler uses up-to-date versions of its source databases.
#             The exact versions change over time. Check g:Profiler website for current data versions.
#             gprofiler2 package version:
gprofiler2_version <- packageVersion("gprofiler2")
cat(paste0("\nQ_ORA_2b: g:Profiler2 R package version used: ", gprofiler2_version, ".\n"))
cat("           Annotation data versions are managed by the g:Profiler service and are regularly updated.\n")

# Q_ORA_3a: How many genesets were returned with what thresholds?
#             Thresholds for DE genes fed into ORA: adj.P.Val < 0.05, |logFC| > 1.0.
#             Threshold for ORA term significance: g:Profiler default (g_SCS corrected p-value < 0.05).
cat("\nQ_ORA_3a: Thresholds for gene selection for ORA: adj.P.Val < ", ORA_ADJP_THRESHOLD, ", |logFC| > ", ORA_LFC_THRESHOLD, ".\n")
cat("           Threshold for ORA term significance: g:Profiler's g_SCS method, effective p-value < 0.05.\n")

cat("           Number of significant gene sets returned for 'all DE genes': ")
if (exists("gostres_all_de") && !is.null(gostres_all_de)) cat(nrow(gostres_all_de$result), "\n") else cat("Not run or no results.\n")
cat("           Number of significant gene sets returned for 'upregulated DE genes': ")
if (exists("gostres_up_de") && !is.null(gostres_up_de)) cat(nrow(gostres_up_de$result), "\n") else cat("Not run or no results.\n")
cat("           Number of significant gene sets returned for 'downregulated DE genes': ")
if (exists("gostres_down_de") && !is.null(gostres_down_de)) cat(nrow(gostres_down_de$result), "\n") else cat("Not run or no results.\n")


# Q_ORA_4: Comparison of results (all DE vs. up/down separately)
# This is more for interpretation, but we can state the number of genes used in each.
cat("\nQ_ORA_4: Comparison of gene lists used for ORA:\n")
cat(paste0("  - All DE genes list size: ", length(all_de_genes_ora_hgnc), "\n"))
cat(paste0("  - Upregulated DE genes list size: ", length(up_regulated_ora_hgnc), "\n"))
cat(paste0("  - Downregulated DE genes list size: ", length(down_regulated_ora_hgnc), "\n"))
cat("  Analyzing up/down lists separately can reveal distinct biological processes affected by up- vs. down-regulation.\n")
cat("  The 'all DE genes' list might highlight overarching themes or pathways where both up/down regulation contribute.\n")

cat("--- COMPLETED SECTION 5 ---\n\n")


# ---
# SECTION 6: INTERPRETATION AND DISCUSSION (Placeholder for Rmd text)
# ---
cat("--- STARTING SECTION 6: INTERPRETATION AND DISCUSSION (Context) ---\n")
cat("This section in the Rmd contained textual interpretation of results.\n")
cat("Key points from Rmd's interpretation:\n")
cat("  - ORA results were compared to the original paper (Oh, Hyunjung et al., 2022).\n")
cat("  - Number of DE genes found here vs. paper: ~608 vs. 835 total; similar ballpark for up/down.\n")
cat("  - Differences attributed to filtering methods and tools (e.g., ClueGo in paper).\n")
cat("  - Evidence from other publications supporting link between genes, depression, and related pathways (Somatostatin, glial cell density, inflammation/apoptosis).\n")
cat("--- COMPLETED SECTION 6 ---\n\n")


# ---
# SECTION 7: CITATIONS (Programmatic parts)
# ---
cat("--- STARTING SECTION 7: CITATIONS ---\n")
cat("Printing citation information for key packages used:\n\n")

# Helper function to print citation if package is loaded
print_citation_if_loaded <- function(pkg_name) {
  if (pkg_name %in% .packages()) {
    cat("Citation for package:", pkg_name, "\n")
    print(citation(pkg_name))
    cat("\n---\n\n")
  } else {
    cat("Package", pkg_name, "not loaded, skipping citation.\n\n")
  }
}

print_citation_if_loaded("tidyverse")
print_citation_if_loaded("edgeR")
# print_citation_if_loaded("GEOmetadb") # Not directly called in final script, but was in RMD setup
print_citation_if_loaded("RColorBrewer")
print_citation_if_loaded("ggplot2")
# print_citation_if_loaded("readxl") # Not directly called in final script
print_citation_if_loaded("dplyr")
print_citation_if_loaded("AnnotationDbi")
print_citation_if_loaded("limma")
print_citation_if_loaded("Biobase")
# print_citation_if_loaded("BiocManager") # Tool, not analysis package
print_citation_if_loaded("biomaRt")
print_citation_if_loaded("magrittr")
print_citation_if_loaded("GEOquery")
# print_citation_if_loaded("RSQLite") # Used by GEOmetadb often
# print_citation_if_loaded("org.Hs.eg.db") # Data package
# print_citation_if_loaded("vegan") # Not directly used in this DGE/ORA workflow
print_citation_if_loaded("gprofiler2")
print_citation_if_loaded("ComplexHeatmap")
print_citation_if_loaded("ggrepel")

cat("--- COMPLETED SECTION 7 ---\n\n")

cat("====== R SCRIPT EXECUTION FINISHED ======\n")
