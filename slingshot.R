#applying slingshot to each seurat object 

#load necessary packages
library("slingshot")
library("Seurat")
library("profmem")
library("tradeSeq")

#specify directory where input data is stored
data_dir = "/home/estagaman/benchmarking_project/data/trajectory/with_DE_genes/seurat/"

#specify desired output directory
out_dir = "/home/estagaman/benchmarking_project/test_simulated_data/slingshot_umap/"

#specify the unique identifier within each input file name 
#in this case, my simulated files are formatted like so: "seurat_glob1.rds", "seurat_glob2.rds", "seurat_1.rds"
#so for each file, the unique identifier would be between the underscore and dot
file_list <- c("glob1", "glob2", "glob3", "glob4", "glob5", "glob6", "1", "2", "3", "4", "5", "6")

#set random seed for reproducibility
set.seed(123)

#function to perform pseudotime ordering for each file
pseudotime_ordering <- function(seurat_object, file, out_dir){

    #save the clustering and embedding information for lineage tracing 
    dimred <- seurat_object@reductions$umap@cell.embeddings #extract umap embedding
    clustering <- seurat_object$RNA_snn_res.0.5 #extract clustering information
    counts <- as.matrix(seurat_object@assays$RNA$counts) #extract counts

    #how do we find out which cluster to put as the root?
    #identify the seurat cluster that matches to the beginning of our ground truth pseudotime
    just_first_10 <- subset(seurat_object@meta.data, pseudotime < 10) #first 10% of pseudotime
    clusters_counted <- table(just_first_10$seurat_clusters) #which seurat clusters contain these cells

    max_frequency_index <- which.max(clusters_counted) #find the cluster containing the most of these cells
    most_common_cluster <- names(clusters_counted[max_frequency_index]) #that is the most common cluster

    root_chosen = most_common_cluster #make that our chosen root

    #run cell ordering now that I know the root state

    #start timer
    ordering_start <- Sys.time()

    p <- profmem({ #use profmem to profile memory for this step

        #arrange the cells into lineages
        lineages <- getLineages(data = dimred,
            clusterLabels = clustering, #use the seurat clustering
            start.clus = root_chosen) #define where to start the trajectory

        #do principal curve smoothing with default parameters from tutorial: https://nbisweden.github.io/workshop-archive/workshop-scRNAseq/2020-01-27/labs/compiled/slingshot/slingshot.html
        curves <- getCurves(lineages, thresh = 0.01, stretch = 0.8, allow.breaks = FALSE, shrink = 0.99)

        #compute pseudotime values
        pseudotime <- slingPseudotime(curves)
    })

    #stop the timer
    ordering_end <- Sys.time()

    #calculate memory usage and time
    total_time <- ordering_end - ordering_start
    total_memory <- sum(as.numeric(p$bytes), na.rm = TRUE)

    #create color palette for plotting tree
    pal <- rainbow(length(unique(clustering)))

    #plot the trajectory and save it as png to output folder
    png(paste0(out_dir, "lineage", file, ".png"), width = 800, height = 600, units = "px", res = 100)
    plot(dimred, col = pal[clustering], asp = 1, pch = 16)
    dev.off()

    #ok, now we're testing for differential expression across pseudotime 
    BiocParallel::register(BiocParallel::SerialParam())

    #fit the generalized additive model
    sce <- fitGAM(counts = as.matrix(counts), sds = curves)

    #save the pseudotime assignments
    write.csv(pseudotime, paste0(out_dir, "pseudotime", file, ".csv"))

    stats_df <- data.frame("ordering_time" = total_time, "memory" = total_memory)
    write.csv(stats_df, paste0(out_dir, "stats_", file, ".csv"))
}

#actually run pseudotime on each of our files 
for (i in file_list){
    seurat_object <- readRDS(paste0(data_dir, "seurat_", i, ".rds")) #if your filenames are formatted differently, change the "seurat_" section to match your filename format

    #perform the pseudotime ordering
    pseudotime_ordering(seurat_object, i, out_dir)
}

#done, yay!
