PAGA_folder <- "/home/estagaman/benchmarking_project/test_simulated_data/PAGA_out_umap_30/seurat"
PAGA_prefix = "pseudotime"

#specify output folder
out_dir <- "/home/estagaman/benchmarking_project/test_simulated_data/correlation_results"

original_object <- readRDS("/home/estagaman/benchmarking_project/data/trajectory/with_DE_genes/seurat/seurat_1.rds")

#use the metadata to find true pseudotime 
true_pseudotime <- original_object@meta.data$pseudotime

file = "glob4"

res <- read.csv(paste0(PAGA_folder, "/", PAGA_prefix, file, ".csv"))

#get the inferred pseudotime
inferred_pseudotime <- res[["dpt_pseudotime"]]

#calculate the pearson correlation
cor_coef <- cor(true_pseudotime, inferred_pseudotime, method = "pearson")

#add this to our list of coefficients
cor_values <- c(cor_values, cor_coef)

#add these to a dataframe so we can plot them in a scatterplot
cor_df <- data.frame(true = true_pseudotime, inferred = inferred_pseudotime)

