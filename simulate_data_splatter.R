#simulate data using splatter 

#load all packages
library(splatter)
library(dplyr)
library(tidyverse)
library(Seurat)
library(patchwork)
library(monocle)
library(reshape)
library(hdf5r)
library(SeuratDisk)

# Create baseline Splat parameters
params <- newSplatParams()

#set the parameters for every simulation
params <- setParams(params, 
                    nGenes = 1000, #1000 genes
                    batchCells = 1000, #1000 cells
                    de.prob = 0.2, #probability of differentially expressed gene is 0.2
                    seed = 100 #provide random seed for reproducibility
)

#simulate with no added dropout, this includes pseudotime steps in the trajectory by default
sim_no_dropout <- splatSimulate(params,
                                method = "paths", #we want a trajectory
                                dropout.type = "none", #no dropout added
                                verbose = FALSE)

#function to add "cell-type-specific" dropouts to this simulated data
add_dropouts <- function(sce, mid){ #takes in single cell experiment and dropout midpoints desired
  
  #now, are dividing our cells up into 4 groups by pseudotime
  start_points <- c(0, 25, 50, 75)
  end_points <- c(25, 50, 75, 100)
  
  #initialize dataframe to keep dropout counts
  rows <- rownames(assay(sce, "counts"))
  cols <- colnames(assay(sce, "counts"))
  
  # initialize empty numeric matrix with 1000 genes x 1000 cells
  counts_w_dropouts <- matrix(NA, nrow = length(rows), ncol = length(cols),
                              dimnames = list(rows, cols))
  
  #for each group of cells:
  for (i in 1:4){
    
    #find the cells within one fourth of pseudotime
    group_cells <- sce$Step <= end_points[i] & sce$Step > start_points[i]
    
    #subset single-cell experiment down to just one group
    one_group <- sce[, sce$Step <= end_points[i] & sce$Step > start_points[i]]
    
    #add dropouts with desired midpoint setting
    data.drop <- splatter:::splatSimDropout(one_group, setParam(params, "dropout.mid", mid[i]))
    
    #extract the dropped counts
    dropped_counts <- assay(data.drop, "counts")
    
    #add them to a new dataframe of counts
    counts_w_dropouts[, group_cells] <- dropped_counts
  }
  
  #save these new counts with dropouts to our single-cell experiment object
  assay(sce, "counts") <- counts_w_dropouts

  #return the final single-cell experiment
  return(sce)
}

#function to calculation the proportion of dataset that equals zero
calculate_zeroes <- function(sce){ #takes in single cell experiment object

  #extract the counts
  counts_mat <- assay(sce, "counts")

  #calculate sum of counts equal to zero, then divide by the total number of counts
  percent_zeroes <- sum(counts_mat == 0) / (nrow(counts_mat) * ncol(counts_mat)) * 100

  #return back percent of zeroes in the dataset
  return(percent_zeroes)
}

#set dropout type as "experiment", which is technically global, but we're doing it in groups
params <- setParams(params, dropout.type = "experiment")

#calculate zeroes in original dataset with no dropouts
calculate_zeroes(sim_no_dropout) #14

##add cell-type specific dropouts with different midpoints per group

data.drop_1 <- add_dropouts(sim_no_dropout, mid = c(0.1, 0.15, 0.2, 0.25))
calculate_zeroes(data.drop_1) #25

data.drop_2 <- add_dropouts(sim_no_dropout, mid = c(0.5, 1, 1.25, 1.5))
calculate_zeroes(data.drop_2) #33

data.drop_3 <- add_dropouts(sim_no_dropout, mid = c(1, 1.25, 1.5, 1.75))
calculate_zeroes(data.drop_3) #37

data.drop_4 <- add_dropouts(sim_no_dropout, mid = c(1.5, 1.75, 2, 2.25))
calculate_zeroes(data.drop_4) #43

data.drop_5 <- add_dropouts(sim_no_dropout, mid = c(0.1, 0.3, 0.5, 0.7))
calculate_zeroes(data.drop_5) #27

#add global dropouts at varied rates across the whole dataset- this is much simpler with Splatter
data.drop_glob1 <- splatter:::splatSimDropout(sim_no_dropout, setParam(params, "dropout.mid", 0.25))

data.drop_glob2 <- splatter:::splatSimDropout(sim_no_dropout, setParam(params, "dropout.mid", 0.5))

data.drop_glob3 <- splatter:::splatSimDropout(sim_no_dropout, setParam(params, "dropout.mid", 0.75))

data.drop_glob4 <- splatter:::splatSimDropout(sim_no_dropout, setParam(params, "dropout.mid", 1))

data.drop_glob5 <- splatter:::splatSimDropout(sim_no_dropout, setParam(params, "dropout.mid", 1.25))

data.drop_glob6 <- splatter:::splatSimDropout(sim_no_dropout, setParam(params, "dropout.mid", 1.5))

#create metadata
metadata_table <- data.frame("cell_id" = colnames(assay(sim_no_dropout, "counts")))
metadata_table$pseudotime <- sim_no_dropout$Step

metadata_table <- metadata_table %>% mutate(
  group = case_when(
    pseudotime <= 25 ~ 1,
    pseudotime > 25 & pseudotime <= 50 ~ 2,
    pseudotime > 50 & pseudotime <= 75 ~ 3,
    pseudotime > 75 & pseudotime <= 100 ~ 4))

rownames(metadata_table) <- metadata_table$cell_id

#save whether a gene is DE
de_genes <- rownames(sim_no_dropout)[rowData(sim_no_dropout)$DEFacPath1 != 1]

#create feature info 
feature_table <- data.frame("feature_id" = rownames(assay(sim_no_dropout, "counts")))
feature_table$gene_mean <- rowData(sim_no_dropout)$GeneMean
feature_table <- feature_table %>% mutate(DE = ifelse(feature_id %in% de_genes, TRUE, FALSE))

#ok, now we make a seurat object for each of these simulated counts tables
seurat_g1 <- CreateSeuratObject(counts = assay(sim_no_dropout, "counts"),
                               meta.data = metadata_table, 
                               min.cells = 0, 
                               min.features = 0,
                               project = "benchmarking_trajectory")

seurat_g2 <- CreateSeuratObject(counts = assay(data.drop_1, "counts"),
                               meta.data = metadata_table, 
                               min.cells = 0, 
                               min.features = 0,
                               project = "benchmarking_trajectory")

seurat_g3 <- CreateSeuratObject(counts = assay(data.drop_2, "counts"),
                               meta.data = metadata_table, 
                               min.cells = 0, 
                               min.features = 0,
                               project = "benchmarking_trajectory")

seurat_g4 <- CreateSeuratObject(counts = assay(data.drop_3, "counts"),
                               meta.data = metadata_table, 
                               min.cells = 0, 
                               min.features = 0,
                               project = "benchmarking_trajectory")

seurat_g5 <- CreateSeuratObject(counts = assay(data.drop_4, "counts"),
                               meta.data = metadata_table, 
                               min.cells = 0, 
                               min.features = 0,
                               project = "benchmarking_trajectory")

seurat_g6 <- CreateSeuratObject(counts = assay(data.drop_5, "counts"),
                               meta.data = metadata_table, 
                               min.cells = 0, 
                               min.features = 0,
                               project = "benchmarking_trajectory")

seurat_glob1 <- CreateSeuratObject(counts = assay(data.drop_glob1, "counts"),
                                meta.data = metadata_table, 
                                min.cells = 0, 
                                min.features = 0,
                                project = "benchmarking_trajectory")

seurat_glob2 <- CreateSeuratObject(counts = assay(data.drop_glob2, "counts"),
                                   meta.data = metadata_table, 
                                   min.cells = 0, 
                                   min.features = 0,
                                   project = "benchmarking_trajectory")

seurat_glob3 <- CreateSeuratObject(counts = assay(data.drop_glob3, "counts"),
                                   meta.data = metadata_table, 
                                   min.cells = 0, 
                                   min.features = 0,
                                   project = "benchmarking_trajectory")

seurat_glob4 <- CreateSeuratObject(counts = assay(data.drop_glob4, "counts"),
                                   meta.data = metadata_table, 
                                   min.cells = 0, 
                                   min.features = 0,
                                   project = "benchmarking_trajectory")

seurat_glob5 <- CreateSeuratObject(counts = assay(data.drop_glob5, "counts"),
                                   meta.data = metadata_table, 
                                   min.cells = 0, 
                                   min.features = 0,
                                   project = "benchmarking_trajectory")

seurat_glob6 <- CreateSeuratObject(counts = assay(data.drop_glob6, "counts"),
                                   meta.data = metadata_table, 
                                   min.cells = 0, 
                                   min.features = 0,
                                   project = "benchmarking_trajectory")


#using seurat QC tutorial: https://satijalab.org/seurat/articles/pbmc3k_tutorial.html

seurat_QC <- function(object){
  
  #do some quality control and feature selection
  #VlnPlot(object, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
  
  #plot scatter of these stats against each other
  #FeatureScatter(object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  
  #normalize the data
  object <- NormalizeData(object, normalization.method = "LogNormalize", scale.factor = 10000)
  
  #note: normalized values are stored in pbmc[["RNA"]]$data
  
  
  #now, we are looking for variable features to base our clustering and trajectory on
  object <- FindVariableFeatures(object, selection.method = "vst", nfeatures = 500)
  
  top10 <- head(VariableFeatures(object), 10)
  
  plot1 <- VariableFeaturePlot(object)
  plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
  #plot1 + plot2
  
  #apply scaling
  all.genes <- rownames(object)
  object <- ScaleData(object, features = all.genes)
  
  #perform dimension reduction 
  object <- RunPCA(object, features = VariableFeatures(object = object)) #results stored in seurat_1[["pca"]]
  
  object <- RunICA(object, features = VariableFeatures(object = object)) #results stored in seurat_1[["ica"]]
  
  object <- RunUMAP(object, features = VariableFeatures(object = object)) #results stored in seurat_1[["umap"]]
  
  #ElbowPlot(object) #we really only need 6-7 PCs, let's move forward with 10 
  
  #do nearest neighbors and clustering
  object <- FindNeighbors(object, dims = 1:10)
  object <- FindClusters(object, resolution = 0.5)
  
  return(object)
}

#ok, so what we want to do now is save this as: 
#1) a cell data set, for monocle
#2) an AnnData hda5 file for PAGA
#3) each single cell experiment saved as an rds file

save_CDS <- function(object, name, out_folder){
  
  #convert to cell data set
  cds <- as.CellDataSet(object, reduction = "ica")
  
  #save as RDS file
  saveRDS(cds, paste0(out_folder, "/", name))
  
  return(cds)
}

save_seurat_object <- function(object, name, out_folder){
  saveRDS(object, file = paste0(out_folder, "/", name))
}

#perform seurat QC
seurat_g1 <- seurat_QC(seurat_g1)
seurat_g2 <- seurat_QC(seurat_g2)
seurat_g3 <- seurat_QC(seurat_g3)
seurat_g4 <- seurat_QC(seurat_g4)
seurat_g5 <- seurat_QC(seurat_g5)
seurat_g6 <- seurat_QC(seurat_g6)
seurat_glob1 <- seurat_QC(seurat_glob1)
seurat_glob2 <- seurat_QC(seurat_glob2)
seurat_glob3 <- seurat_QC(seurat_glob3)
seurat_glob4 <- seurat_QC(seurat_glob4)
seurat_glob5 <- seurat_QC(seurat_glob5)
seurat_glob6 <- seurat_QC(seurat_glob6)

#convert to CDS for Monocle
cds_folder <- "/Users/elise/Downloads/BIOI_benchmarking_project/with_DE_genes/CDS"
cds_1 <- save_CDS(seurat_g1, "cds1.rds", cds_folder)
cds_2 <- save_CDS(seurat_g2, "cds2.rds", cds_folder)
cds_3 <- save_CDS(seurat_g3, "cds3.rds", cds_folder)
cds_4 <- save_CDS(seurat_g4, "cds4.rds", cds_folder)
cds_5 <- save_CDS(seurat_g5, "cds5.rds", cds_folder)
cds_6 <- save_CDS(seurat_g6, "cds6.rds", cds_folder)

cds_glob1 <- save_CDS(seurat_glob1, "cds_glob1.rds", cds_folder)
cds_glob2 <- save_CDS(seurat_glob2, "cds_glob2.rds", cds_folder)
cds_glob3 <- save_CDS(seurat_glob3, "cds_glob3.rds", cds_folder)
cds_glob4 <- save_CDS(seurat_glob4, "cds_glob4.rds", cds_folder)
cds_glob5 <- save_CDS(seurat_glob5, "cds_glob5.rds", cds_folder)
cds_glob6 <- save_CDS(seurat_glob6, "cds_glob6.rds", cds_folder)

seurat_folder <- "/Users/elise/Downloads/BIOI_benchmarking_project/with_DE_genes/seurat"
save_seurat_object(seurat_g1, "seurat_1.rds", seurat_folder)
save_seurat_object(seurat_g2, "seurat_2.rds", seurat_folder)
save_seurat_object(seurat_g3, "seurat_3.rds", seurat_folder)
save_seurat_object(seurat_g4, "seurat_4.rds", seurat_folder)
save_seurat_object(seurat_g5, "seurat_5.rds", seurat_folder)
save_seurat_object(seurat_g6, "seurat_6.rds", seurat_folder)

save_seurat_object(seurat_glob1, "seurat_glob1.rds", seurat_folder)
save_seurat_object(seurat_glob2, "seurat_glob2.rds", seurat_folder)
save_seurat_object(seurat_glob3, "seurat_glob3.rds", seurat_folder)
save_seurat_object(seurat_glob4, "seurat_glob4.rds", seurat_folder)
save_seurat_object(seurat_glob5, "seurat_glob5.rds", seurat_folder)
save_seurat_object(seurat_glob6, "seurat_glob6.rds", seurat_folder)

save_h5ad <- function(object, name, out_folder){
  
  #follow instructions for converting seurat object from v5 to v3, then save using sceasy
  #https://github.com/satijalab/seurat/issues/8220
  
  #convert from v5 to v3 assay
  object[["RNA3"]] <- as(object = object[["RNA"]], Class = "Assay")
  DefaultAssay(object) <- "RNA3"
  object[["RNA"]] <- NULL
  object <- RenameAssays(object = object, RNA3 = 'RNA')
  
  #make the path we want to save this file to
  path_to_save <- paste0(out_folder, "/", name)
  
  #save as SeuratH5
  SaveH5Seurat(object, filename = path_to_save)
  
  #convert to h5ad
  Convert(path_to_save, dest = "h5ad", overwrite = "T")
}

h5ad_out_folder <- "/Users/elise/Downloads/BIOI_benchmarking_project/with_DE_genes/h5ad"
save_h5ad(seurat_g1, "seurat_1.h5Seurat", h5ad_out_folder)
save_h5ad(seurat_g2, "seurat_2.h5Seurat", h5ad_out_folder)
save_h5ad(seurat_g3, "seurat_3.h5Seurat", h5ad_out_folder)
save_h5ad(seurat_g4, "seurat_4.h5Seurat", h5ad_out_folder)
save_h5ad(seurat_g5, "seurat_5.h5Seurat", h5ad_out_folder)
save_h5ad(seurat_g6, "seurat_6.h5Seurat", h5ad_out_folder)

save_h5ad(seurat_glob1, "seurat_glob1.h5Seurat", h5ad_out_folder)
save_h5ad(seurat_glob2, "seurat_glob2.h5Seurat", h5ad_out_folder)
save_h5ad(seurat_glob3, "seurat_glob3.h5Seurat", h5ad_out_folder)
save_h5ad(seurat_glob4, "seurat_glob4.h5Seurat", h5ad_out_folder)
save_h5ad(seurat_glob5, "seurat_glob5.h5Seurat", h5ad_out_folder)
save_h5ad(seurat_glob6, "seurat_glob6.h5Seurat", h5ad_out_folder)

#save R workspace
save.image(file = "/Users/elise/Downloads/BIOI_benchmarking_project/with_DE_genes/simulation_workspace.RData")

#save features information with DE
write.csv(feature_table, "/Users/elise/Downloads/BIOI_benchmarking_project/with_DE_genes/feature_info.csv")

