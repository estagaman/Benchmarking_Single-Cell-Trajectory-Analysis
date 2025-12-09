#measure the correlation between actual pseudotime and assigned pseudotime by the tool 

library(patchwork)
library(ggplot2)

#which tool would we like to assess? put its output folder here

#Monocle
monocle_folder <- "/home/estagaman/benchmarking_project/test_simulated_data/monocle_DDR_mem"
monocle_prefix = "results_"

#PAGA
PAGA_folder <- "/home/estagaman/benchmarking_project/test_simulated_data/PAGA_out_umap_30/seurat"
PAGA_prefix = "pseudotime_"

#Slingshot
Slingshot_folder <- "/home/estagaman/benchmarking_project/test_simulated_data/slingshot_out_umap"
Slingshot_prefix = "pseudotime"

#do we want to look at global dropouts or cell-type-specific dropouts?
check <- "cell_type"

#specify output folder
out_dir <- "/home/estagaman/benchmarking_project/test_simulated_data/correlation_results"

if (check == "global"){

    file_list <- c("1", "glob1", "glob2", "glob3", "glob4", "glob5", "glob6")
    labels <- c("no dropout", "0.25 dropout", "0.5 dropout", "0.75 dropout", "1 dropout", "1.25 dropout", "1.5 dropout")

} else if (check == "cell_type"){

    file_list <- c("1", "2", "3", "4", "5", "6")
    labels <- c("no dropout", "[0.1, 0.15, 0.2, 0.25]", "[0.5, 1, 1.25, 1.5]", "[1, 1.25, 1.5, 1.75]", "[1.5, 1.75, 2, 2.25]", "[0.1, 0.3, 0.5, 0.7]")

}

#step 1: load in the actual pseudotime information

#we're going to do this with the Seurat object 
original_object <- readRDS("/home/estagaman/benchmarking_project/data/trajectory/with_DE_genes/seurat/seurat_1.rds")

#use the metadata to find true pseudotime 
true_pseudotime <- original_object@meta.data$pseudotime

#step 2, for each file, compute pearson correlation coefficient

find_cor <- function(file_list, prefix, results_folder, col_name){

    cor_values <- c()
    plot_list <- list()

    for (file in file_list){

        #read in the results file
        res <- read.csv(paste0(results_folder, "/", prefix, file, ".csv"))

        #get the inferred pseudotime
        inferred_pseudotime <- res[[col_name]]

        #calculate the pearson correlation
        cor_coef <- cor(true_pseudotime, inferred_pseudotime, method = "pearson")

        #add this to our list of coefficients
        cor_values <- c(cor_values, cor_coef)

        #add these to a dataframe so we can plot them in a scatterplot
        cor_df <- data.frame(true = true_pseudotime, inferred = inferred_pseudotime)

        #create a scatterplot
        plot <- ggplot(data = cor_df, aes(x = true, y = inferred)) + 
                geom_point(color = "darkblue") +
                annotate("text",
                    x = Inf, y = -Inf,         # bottom-right corner
                    label = paste0("r = ", round(cor_coef, 3)),
                    hjust = 1.1, vjust = -0.5, # adjust inward
                    size = 7
                ) + 
                theme(aspect.ratio = 1)
    
        #add this to our list of plots
        plot_list <- append(plot_list, list(plot))
    }

    return(list(plot_list, cor_values))

}

#initialize correlation df
correlation_df <- data.frame(placeholder = c(rep(NA, length(file_list))), row.names = labels)

#run results for monocle
monocle_results <- find_cor(file_list, monocle_prefix, monocle_folder, "Pseudotime")
correlation_df$monocle = monocle_results[[2]]

#now we can make a plot with all the scatterplots combined together
monocle_scatter <- wrap_plots(monocle_results[[1]], ncol = 6)

#save monocle scatterplots
ggsave(paste0(out_dir, "/monocle_DDR_scatter.png"), monocle_scatter, height = 5, width = 30)


#run results for PAGA
PAGA_results <- find_cor(file_list, PAGA_prefix, PAGA_folder, "dpt_pseudotime")
correlation_df$PAGA = PAGA_results[[2]]

#now we can make a plot with all the scatterplots combined together
paga_scatter <- wrap_plots(PAGA_results[[1]], ncol = 6)

#save monocle scatterplots
ggsave(paste0(out_dir, "/PAGA_umap_scatter.png"), paga_scatter, height = 5, width = 30)

#run results for Slingshot
SS_results <- find_cor(file_list, Slingshot_prefix, Slingshot_folder, "Lineage1")
correlation_df$Slingshot = SS_results[[2]]

#now we can make a plot with all the scatterplots combined together
ss_scatter <- wrap_plots(SS_results[[1]], ncol = 6)

#save monocle scatterplots
ggsave(paste0(out_dir, "/Slingshot_umap_scatter.png"), ss_scatter, height = 5, width = 30)

correlation_df$placeholder = NULL
write.csv(correlation_df, paste0(out_dir, "/cor_values_", "global_default.csv"))

final_plot <- monocle_scatter / paga_scatter / ss_scatter

ggsave(paste0(out_dir, "/cell_type_default_scatter.png"), final_plot, height = 15, width = 30)
