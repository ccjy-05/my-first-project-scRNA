# Single-cell Transcriptomic Analysis of Spermatogenesis and Virtual Knockout

This repository contains the R scripts for single-cell RNA sequencing (scRNA-seq) analysis of mouse testis tissue. The pipeline includes data quality control, clustering, cell type annotation across the spermatogenesis trajectory, and downstream virtual knockout simulations (using `scTenifoldKnk`) for target genes (e.g., *stmn1a* and *tubb2b*).

##  Data Source & Library Preparation

**Biological Sample & Library Prep**:
Single-cell RNA sequencing libraries were prepared using the DNBelab C Series High-throughput Single-Cell RNA Library Preparation Kit V2.0 (MGI, China) according to the manufacturer's instructions. Briefly, microfluidic droplet generation was performed on the DNBelab C4 platform, followed by reverse transcription, cDNA amplification, and library construction. High-throughput sequencing of the libraries was carried out by Shenzhen Geno-Immune Medical Institute.

**Data Processing**:
Raw sequencing data (FASTQ) were processed using the **DNBC4tools (v2.1)**, the official analysis pipeline for the DNBelab C Series platform. 
*Note*: The output expression matrix from DNBC4tools was formatted into a 10x Genomics-compatible structure (`barcodes.tsv.gz`, `features.tsv.gz`, `matrix.mtx.gz`), allowing for direct loading into R via `Seurat::Read10X()`.

##  Dependencies
The code is written in R (>= 4.0). Please install the required packages before running the scripts:

```R
# Install standard CRAN / Bioconductor packages
install.packages(c("Seurat", "ggplot2", "dplyr", "tidyr", "patchwork", "Matrix", "parallel"))
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("rhdf5"))

# Install scTenifoldNet and scTenifoldKnk from GitHub
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("cailab-tamu/scTenifoldNet")
devtools::install_github("cailab-tamu/scTenifoldKnk")