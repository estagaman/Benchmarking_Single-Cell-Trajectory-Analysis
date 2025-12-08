# Benchmarking Tools for Single-Cell Trajectory Analysis

## Introduction: 

Today, single-cell RNA-sequencing allows for profiling of gene expression in thousands of cells simultaneously. Single-cell trajectory analysis tools aim to situate these cells within a timecourse, inferring cell differentiation and evolutionary processes. However, data at single-cell resolution presents unique obstacles, particularly low counts, dropouts, and cell-type-specific biological effects. In this study, I test the efficiency of common trajectory analysis tools, including Monocle2, Slingshot, and PAGA, as well as their robustness to both global and cell-type-specific dropouts. Analysis reveals that both PAGA and Slingshot most accurately infer true pseudotime in the presence of dropouts, with superior time and memory efficiency compared to Monocle2. This implies that improvements upon the original minimum-spanning-tree based algorithms used by Monocle, including principal curve smoothing, clustering-based approaches, and connectivity graphing, notably improve stability and efficiency. Future single-cell analysis tools may leverage this knowledge to develop new graph-based tools. 

## How to Simulate Data: 

All data simulation code for this project is included in script data_simulation.R

Necessary R packages include: 
  - Splatter
  - Seurat
  - Monocle
  - SeuratDisk
  - tidyverse
  - dplyr

## Running Each Tool 

### Tool #1: Monocle2 

Code to run this tool is available in script Monocle2.R

Input: 
  - CellDataSet object including raw counts from scRNA-seq, saved as .rds file
  - vector prefixes for each file (section that comes before ".rds" (if using my test data, keep prefixes the same
  - desired output folder name

Output: 
  - Stats: time and memory usage statistics for each file
  - 


