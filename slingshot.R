#applying slingshot to each seurat object 

#load slingshot 
library("slingshot")
library("Seurat")
library("profmem")
library("tradeSeq")

data_dir = "/home/estagaman/benchmarking_project/data/trajectory/with_DE_genes/seurat/"
out_dir = "/home/estagaman/benchmarking_project/test_simulated_data/slingshot_out_ica/"

file_list <- c("glob1", "glob2", "glob3", "glob4", "glob5", "glob6", "1", "2", "3", "4", "5", "6")

set.seed(123)

pseudotime_ordering <- function(seurat_object, file, out_dir){

    #save the clustering and embedding information for lineage tracing 
    dimred <- seurat_object@reductions$ica@cell.embeddings
    clustering <- seurat_object$RNA_snn_res.0.5
    counts <- as.matrix(seurat_object@assays$RNA$counts)

    set.seed(123)

    #how do we find out which cluster to put as the root?
    #identify the seurat cluster that matches to the beginning of our ground truth pseudotime
    just_first_10 <- subset(seurat_object@meta.data, pseudotime < 10)
    clusters_counted <- table(just_first_10$seurat_clusters)

    max_frequency_index <- which.max(clusters_counted)
    most_common_cluster <- names(clusters_counted[max_frequency_index])

    root_chosen = most_common_cluster

    #run it again now that I know where I want it to end 
    ordering_start <- Sys.time()

    p <- profmem({ 
        lineages <- getLineages(data = dimred,
            clusterLabels = clustering,
            start.clus = root_chosen) #define where to start the trajectories

        curves <- getCurves(lineages, thresh = 0.01, stretch = 0.8, allow.breaks = FALSE, shrink = 0.99)

        pseudotime <- slingPseudotime(curves)
    })

    ordering_end <- Sys.time()

    #calculate memory usage and time
    total_time <- ordering_end - ordering_start
    total_memory <- sum(as.numeric(p$bytes), na.rm = TRUE)

    pal <- rainbow(length(unique(clustering)))

    png(paste0(out_dir, "lineage", file, ".png"), width = 800, height = 600, units = "px", res = 100)
    plot(dimred, col = pal[clustering], asp = 1, pch = 16)
    dev.off()

    #ok, now we're testing for differential expression across pseudotime 
    BiocParallel::register(BiocParallel::SerialParam())

    #fit the generalized additive model
    sce <- fitGAM(counts = as.matrix(counts), sds = curves)

    #test for association with pseudotime
    pseudotime_association <- associationTest(sce)
    pseudotime_association$fdr <- p.adjust(pseudotime_association$pvalue, method = "fdr")
    pseudotime_association <- pseudotime_association[order(pseudotime_association$pvalue), ]
    pseudotime_association$feature_id <- rownames(pseudotime_association)

    #save these results with p-values
    write.csv(pseudotime_association, paste0(out_dir, "pseudotime_DE", file, ".csv"))

    #save the pseudotime assignments
    write.csv(pseudotime, paste0(out_dir, "pseudotime", file, ".csv"))

    stats_df <- data.frame("ordering_time" = total_time, "memory" = total_memory)

    write.csv(stats_df, paste0(out_dir, "stats_", file, ".csv"))

}

#actually run pseudotime on each of our files 
for (i in file_list){
    seurat_object <- readRDS(paste0(data_dir, "seurat_", i, ".rds"))

    pseudotime_ordering(seurat_object, i, out_dir)
}

#done, yay!
