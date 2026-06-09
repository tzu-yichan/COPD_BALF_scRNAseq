# COPD_BALF_scRNAseq

Code repository accompanying:

**Lipid Alterations in Alveolar Surfactant Compromise Alveolar Macrophage Immune Function in COPD**

## Study Overview

This repository contains scripts and metadata used for:

- Single-cell RNA sequencing analysis of BALF cells
- Quality control and data integration
- Cell type annotation using SingleR
- Macrophage re-clustering
- Differential expression analysis
- Hallmark gene set enrichment analysis (GSEA)
- Visualization of single-cell transcriptomic data

## Study Cohorts

- Control (n = 5)
- COPD-quit (n = 4)
- COPD-active (n = 5)

Bronchoalveolar lavage fluid (BALF) cells from subjects within each clinical group were pooled prior to library preparation and sequenced using the 10x Genomics Fixed RNA Profiling platform.

## Repository Structure

```text
scripts/   Analysis scripts
metadata/  Cell-level annotations and metadata
figures/   Figure resources
