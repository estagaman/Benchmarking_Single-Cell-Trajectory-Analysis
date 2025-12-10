# Benchmarking Tools for Single-Cell Trajectory Analysis

## Introduction: 

Today, single-cell RNA-sequencing allows for profiling of gene expression in thousands of cells simultaneously. Single-cell trajectory analysis tools aim to situate these cells within a time course, inferring cell differentiation and evolutionary processes. However, data at single-cell resolution presents unique obstacles, particularly low counts, dropouts, and cell-type-specific biological effects. In this study, I test the efficiency of common trajectory analysis tools, including Monocle2, Slingshot, and PAGA, as well as their robustness to both global and cell-type-specific dropouts. Analysis reveals that both PAGA and Slingshot more accurately infer true pseudotime in the presence of moderate, global dropout rates. These tools also required only 0-4 minutes to execute cell ordering, whereas Monocle2 required more than 20 seconds total. 
However, Monocle2 performed more consistently with severe, cell-type-specific sparsity and showed up to twice as much sensitivity in identifying differentially expressed genes. This implies that improvements upon the original minimum-spanning-tree based algorithms used by Monocle, including principal curve smoothing, clustering-based approaches, and connectivity graphing, improve efficiency, but not necessarily accuracy or sensitivity in the presence of noise. Future single-cell analysis tools may combine Monocle2’s advantages in precision, recall, and noise handling with the computational efficiency of PAGA and Slingshot. Currently, recommended trajectory analysis tools depends on the needs of the particular datasets, amount of dropout heterogeneity, computational resources, and specific research question. 

## Required Software: 

R: version 3.0 or later
    * download link: https://www.r-project.org 

Python: version 3.9 or later
    * download link: https://www.python.org/downloads/

Anaconda/Miniconda: 
    * download link: [https://www.python.org/downloads/](https://docs.conda.io/projects/conda/en/stable/user-guide/install/index.html)

### Required R packages: 
    * Bioconductor
    * tidyverse
    * monocle (version 2.0+)
    * splatter
    * dplyr
    * Seurat
    * patchwork
    * reshape
    * SeuratDisk
    * slingshot
    * tradeSeq
    * profmem
    * igraph
    * ggplot2
    * patchwork
  
### Required Python modules: 
    * pandas
    * numpy
    * matplotlib
    * scanpy (must be activated in conda environment to run PAGA)
    * scipy
    * time
    * os
    * tracemalloc

## How to Simulate Data: 

All data simulation code for this project is included in script data_simulation.R

## Running Each Tool:

### Tool #1: Monocle2 

Code to run this tool is available in script Monocle.R

Files used as input for this script should be in the same folder, data_dir, with filename structure "cds_{file_ID}.rds" where file_ID is a unique identifier for that dataset

Input: 
  - data_dir: folder of CellDataSet objects, one for each dataset, including raw counts from scRNA-seq, saved as .rds file
  - file_list: vector of identifiers (file_IDs) for each file (section of file name that comes before ".rds" (if using my simulated data, keep file_list the same)
  - out_dir: desired output folder path

Output: 
  - stats_{file_ID}.csv: time and memory usage statistics for each dataset
  - pseudotime_{file_ID}.csv: pseudotime assignments for each dataset
  - tree_{file_ID}.csv: plot of trajectory for visual evaluation

### Tool #2: PAGA

Code to run this tool is available in script PAGA.py

Files used as input for this script should be in the same folder, data_dir, with filename structure "seurat_{file_ID}.h5ad" where file_ID is a unique identifier for that dataset

Input: 
  - data_dir: folder of h5ad files, one for each dataset, including raw counts from scRNA-seq and UMAP projection and clustering computed by Seurat
  - file_list: vector of identifiers (file_IDs) for each file (section of file name that comes before ".h5ad" (if using my simulated data, keep file_list the same)
  - out_dir: desired output folder path

Output: 
  - stats_{file_ID}.csv: time and memory usage statistics for each dataset
  - pseudotime_{file_ID}.csv: pseudotime assignments for each dataset

### Tool #3: Slingshot

Code to run this tool is available in script slingshot.R

Files used as input for this script should be in the same folder, data_dir, with filename structure "seurat_{file_ID}.rds" where file_ID is a unique identifier for that dataset

Input: 
  - data_dir: folder of .rds files, one for each dataset, including raw counts from scRNA-seq and UMAP projection and clustering collected in a Seurat object
  - file_list: vector of identifiers (file_IDs) for each file (section of file name that comes before ".rds" (if using my simulated data, keep file_list the same)
  - out_dir: desired output folder path

Output: 
  - stats_{file_ID}.csv: time and memory usage statistics for each dataset
  - pseudotime_{file_ID}.csv: pseudotime assignments for each dataset

## Calculating Differentially Expressed Genes Across Pseudotime

For the purpose of continuity between tools, pseudotime results from all tools were tested for differentially expressed genes using the StartVsEndTest() from tradeSeq R package. Code is available for this in tradeseq_start_end.R 

Results of each tool must be run individually in tradeSeq, with parameters out_dir and time_column adjusted like so: 

```{r}
#specify the directory where your pseudotime results are stored for a particular tool 
out_dir = "/home/estagaman/benchmarking_project/test_simulated_data/monocle_DDR"

#give the prefix assigned to the csv with pseudotime values. If using my code, this should be "pseudotime_"
file_prefix = "/pseudotime_"

#give the column name where inferred pseudotime is stored. 
time_column = "Pseudotime"
    #for monocle: "Pseudotime"
    #for PAGA: "dpt_pseudotime" 
    #for slingshot: "Lineage"
```
Results will output to the out_dir specified, formatted as pseudotimeDE_{file_ID}.csv

## Generating Benchmarking Results 

The following section is ordered according to the figures in my benchmarking paper results: 

### Figure 1: Time and Space Complexity for Global Dropout Datasets: 

This figure was generated using evaluation_code/time_space.R, with the following user settings: 

```{r}

#folder with monocle time and space statistics, with beginning prefix of all filenames with statistics
monocle_folder <- "/home/estagaman/benchmarking_project/test_simulated_data/monocle_DDR"
monocle_prefix = "stats_"

#folder with PAGA time and space statistics, with beginning prefix of all filenames with statistics
PAGA_folder <- "/home/estagaman/benchmarking_project/test_simulated_data/PAGA_umap"
PAGA_prefix = "stats_"

#folder with slingshot time and space statistics, with beginning prefix of all filenames with statistics
Slingshot_folder <- "/home/estagaman/benchmarking_project/test_simulated_data/slingshot_umap"
Slingshot_prefix = "stats_"

#do we want to look at global dropouts or cell-type-specific dropouts?
check <- "global"

#specify output folder
out_dir <- "/home/estagaman/benchmarking_project/test_simulated_data/time_space_results"

```

Folders and output directory should be customized according to where you stored your results from initial runs of each tool. 

The final plot for Figure 1 is called "time_space_global_default.png", located inside the specified output directory. 

### Table 2: Space Complexity for Global Dropout Datasets 

This information is calculated internally with the same script and settings as Figure 1, evaluation_code/time_space.R

The final output file is called "time_space_global_default.csv", located inside the specified output directory. 

The resulting metrics were manually transferred into the paper. 


### Table 3: Memory Usage for Global Dropout Datasets 

This information is contained in the same output file as Table 2, "time_space_global_default.csv", and generated with script evaluation_code/time_space.R. 

Resulting metrics were manually transferred into the paper. 


### Figure 2: Accuracy Measures in Global Dropout Datasets

#### part a) generated using evaluation_code/pt_correlation.R, with the following user settings: 

```{r}
#Monocle-specify folder with results and prefix of files with pseudotime assignments
monocle_folder <- "/home/estagaman/benchmarking_project/test_simulated_data/monocle_DDR"
monocle_prefix = "pseudotime_"

#PAGA-specify folder with results and prefix of files with pseudotime assignments
PAGA_folder <- "/home/estagaman/benchmarking_project/test_simulated_data/PAGA_umap"
PAGA_prefix = "pseudotime_"

#Slingshot-specify folder with results and prefix of files with pseudotime assignments
Slingshot_folder <- "/home/estagaman/benchmarking_project/test_simulated_data/slingshot_umap"
Slingshot_prefix = "pseudotime_"

#do we want to look at global dropouts or cell-type-specific dropouts?
check <- "global"

#specify output folder to save results
out_dir <- "/home/estagaman/benchmarking_project/test_simulated_data/correlation_results"

#give path to seurat object with no dropouts added: 
original_seurat_object_path <- "/home/estagaman/benchmarking_project/data/trajectory/with_DE_genes/seurat/seurat_1.rds"
```

Folders and output directory should be customized according to where you stored your results from initial runs of each tool. Parameter check should be set to "global" for global dropout datasets. 

The final plot for Figure 2a is called "global_default_scatter.png", located inside the specified output directory.

#### part b) generated using evaluation_code/prec_recall.R, with the following user settings: 

```{r}
#folder with monocle results, and beginning of filenames with differential expression results
monocle_folder <- "/home/estagaman/benchmarking_project/test_simulated_data/monocle_DDR"
monocle_prefix = "pseudotime_DE"

#folder with monocle results, and beginning of filenames with differential expression results
PAGA_folder <- "/home/estagaman/benchmarking_project/test_simulated_data/PAGA_umap"
PAGA_prefix = "pseudotime_DE" #using the start vs end test

#folder with monocle results, and beginning of filenames with differential expression results
Slingshot_folder <- "/home/estagaman/benchmarking_project/test_simulated_data/slingshot_umap"
Slingshot_prefix = "pseudotime_DE"

#do we want to look at global dropouts or cell-type-specific dropouts?
check <- "global"

#do we want to use p-value or FDR-corrected p-value as the threshold?
threshold <- "pvalue"

#specify output folder
out_dir <- "/home/estagaman/benchmarking_project/test_simulated_data/prec_recall_results"

#provide a csv with a column of gene IDs and one column on whether the gene is truly differentially expressed
DE_genes <- read.csv("/home/estagaman/benchmarking_project/data/trajectory/with_DE_genes/feature_info.csv")
```

Folders and output directory should be customized according to where you stored your results from initial runs of each tool. Parameter check should be set to "global" for global dropout datasets. 

The final plot for Figure 2b is called "default_prec_recallglobal.png", located inside the specified output directory.

### Figure 3: Accuracy Measures in Cell-Type-Specific Dropout Datasets

#### part a) generated using evaluation_code/pt_correlation.R, with the following changes in settings from Figure 2: 

```{r}
#do we want to look at global dropouts or cell-type-specific dropouts?
check <- "cell_type"

```

Parameter check should be set to "cell_type" for cell-type-specific dropout datasets. 

The final plot for Figure 3a is called "cell_type_default_scatter.png", located inside the specified output directory.

#### part b) generated using evaluation_code/prec_recall.R, with the following changes in settings from Figure 2: 

```{r}
#do we want to look at global dropouts or cell-type-specific dropouts?
check <- "cell_type"

```

Parameter check should be set to "cell_type" for cell-type-specific dropout datasets. 

The final plot for Figure 3b is called "default_prec_recallcell_type.png", located inside the specified output directory.
