#load necessary packages
library("tradeSeq")
library("monocle")

#give directory where original CellDataSets were stored for Monocle specifically
data_dir = "/home/estagaman/benchmarking_project/data/trajectory/with_DE_genes/CDS"

#give output directory where pseudotime assignments were saved for the tool you're testing
#options based on previous code: monocle_DDR, PAGA_umap, or slingshot_umap
out_dir = "/home/estagaman/benchmarking_project/test_simulated_data/monocle_DDR" #in this case I'm running monocle

file_prefix = "/pseudotime_" #give the prefix assigned to the csv with pseudotime values
time_column = "Pseudotime" #give the column name where inferred pseudotime is stored. 
    #for monocle: "Pseudotime:"
    #for PAGA: "dpt_pseudotime" 
    #for slingshot: "Lineage"

#give a list of the unique file identifiers associated with each dataset
#with this data, it will read in filenames with your specified prefix + file from file_list + .csv
#for example: "pseudotime_glob1.csv", "pseudotime_glob2.csv", "pseudotime_1.csv"
file_list <- c("glob1", "glob2", "glob3", "glob4", "glob5", "glob6", "1", "2", "3", "4", "5", "6")

#set a seed for reproducibility
set.seed(123)

#set the parameters
BiocParallel::register(BiocParallel::SerialParam())

#for each file we want to test:
for (i in file_list){

    #read the original CellDataSet for that data
    cds <- readRDS(paste0(data_dir, "/cds_", i, ".rds"))

    #extract the counts
    counts <- exprs(cds)

    #read in the pseudotime assignments from the tool indicated
    pseudotime = read.csv(paste0(out_dir, file_prefix, i, ".csv"))

    #create a dataframe with the Cell ID and the inferred pseudotime
    pseudotime = data.frame(row.names = pseudotime$X, pseudotime = pseudotime[, time_column])

    #create a vector of CellWeights if multiple lineages. In this case, all weights are 1
    cellWeights = data.frame(row.names = rownames(pseudotime), weights = c(rep(1, nrow(pseudotime))))

    #fit the generalized additive model (using 3 knots for smaller sample size
    sce <- fitGAM(counts = as.matrix(counts), pseudotime = pseudotime, cellWeights = cellWeights, nknots = 3)

    #test for association with pseudotime using startVsEndTest
    pseudotime_association <- startVsEndTest(sce)
    pseudotime_association$fdr <- p.adjust(pseudotime_association$pvalue, method = "fdr")
    pseudotime_association <- pseudotime_association[order(pseudotime_association$pvalue), ]
    pseudotime_association$feature_id <- rownames(pseudotime_association)

    #save these results with p-values
    write.csv(pseudotime_association, paste0(out_dir, "/pseudotime_DE", i, ".csv"))

}
