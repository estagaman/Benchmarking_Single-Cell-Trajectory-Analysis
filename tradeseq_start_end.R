library("tradeSeq")
library("monocle")

data_dir = "/home/estagaman/benchmarking_project/data/trajectory/with_DE_genes/CDS"
out_dir = "/home/estagaman/benchmarking_project/test_simulated_data/monocle_DDR_mem"

file_prefix = "/results_"
time_column = "Pseudotime"

file_list <- c("glob1", "glob2", "glob3", "glob4", "glob5", "glob6", "1", "2", "3", "4", "5", "6")

set.seed(123)

BiocParallel::register(BiocParallel::SerialParam())

for (i in file_list){

    cds <- readRDS(paste0(data_dir, "/cds_", i, ".rds"))
    counts <- exprs(cds)

    pseudotime = read.csv(paste0(out_dir, file_prefix, i, ".csv"))
    pseudotime = data.frame(row.names = pseudotime$X, pseudotime = pseudotime[, time_column])

    cellWeights = data.frame(row.names = rownames(pseudotime), weights = c(rep(1, nrow(pseudotime))))

    #fit the generalized additive model
    sce <- fitGAM(counts = as.matrix(counts), pseudotime = pseudotime, cellWeights = cellWeights, nknots = 3)

    #test for association with pseudotime
    pseudotime_association <- startVsEndTest(sce)
    pseudotime_association$fdr <- p.adjust(pseudotime_association$pvalue, method = "fdr")
    pseudotime_association <- pseudotime_association[order(pseudotime_association$pvalue), ]
    pseudotime_association$feature_id <- rownames(pseudotime_association)

    #save these results with p-values
    write.csv(pseudotime_association, paste0(out_dir, "/pseudotime_DEse3", i, ".csv"))

}
