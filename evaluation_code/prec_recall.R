library(patchwork)
library(ggplot2)

#Monocle
monocle_folder <- "/home/estagaman/benchmarking_project/test_simulated_data/monocle_DDR_mem"
#monocle_prefix = "DEG_"
monocle_prefix = "pseudotime_DEse3"

#PAGA
PAGA_folder <- "/home/estagaman/benchmarking_project/test_simulated_data/PAGA_out_umap_30/seurat"
PAGA_prefix = "pseudotime_DEse3" #using the start vs end test

#Slingshot
Slingshot_folder <- "/home/estagaman/benchmarking_project/test_simulated_data/slingshot_out_umap"
Slingshot_prefix = "pseudotime_DEse3"

#do we want to look at global dropouts or cell-type-specific dropouts?
check <- "cell_type"
threshold <- "pvalue"

#specify output folder
out_dir <- "/home/estagaman/benchmarking_project/test_simulated_data/prec_recall_results"

#load in the true DE genes
DE_genes <- read.csv("/home/estagaman/benchmarking_project/data/trajectory/with_DE_genes/feature_info.csv")

if (check == "global"){

    file_list <- c("1", "glob1", "glob2", "glob3", "glob4", "glob5", "glob6")
    labels <- c("no dropout", "0.25 dropout", "0.5 dropout", "0.75 dropout", "1 dropout", "1.25 dropout", "1.5 dropout")

} else if (check == "cell_type"){

    file_list <- c("1", "2", "3", "4", "5", "6")
    labels <- c("no dropout", "[0.1, 0.15, 0.2, 0.25]", "[0.5, 1, 1.25, 1.5]", "[1, 1.25, 1.5, 1.75]", "[1.5, 1.75, 2, 2.25]", "[0.1, 0.3, 0.5, 0.7]")

}

calculate_metrics <- function(file_list, prefix, results_folder, method){

    stats_df <- data.frame(precision = c(), recall = c(), method = c())

    for (file in file_list){

        #read in the results file
        res <- read.csv(paste0(results_folder, "/", prefix, file, ".csv"))

        #make a vector of values TRUE and FALSE that show whether a feature was identified as differentially expressed or not

        if (threshold == "fdr"){
            if ("DE" %in% colnames(res)){
                DE_inferred = res$DE == "DE"
            } else {
                DE_inferred <- ifelse(is.na(res$fdr), FALSE, res$fdr < 0.05)
            }
        } else if (threshold == "pvalue"){ #if we're using pvalue threshold
            if ("pvalue" %in% colnames(res)){
                DE_inferred = res$pvalue < 0.05
            }
        }

        DE_inferred <- factor(DE_inferred, levels = c(FALSE, TRUE))
        DE_truth <- factor(DE_genes$DE, levels = c(FALSE, TRUE))

        #calculate true/false positives/negatives
        outcome_table <- table(inferred = DE_inferred, truth = DE_truth)

        true_pos <- outcome_table["TRUE", "TRUE"]
        false_pos <- outcome_table["TRUE", "FALSE"]
        false_neg <- outcome_table["FALSE", "TRUE"]
        true_neg <- outcome_table["FALSE", "FALSE"]

        precision = true_pos / (true_pos + false_pos)
        recall = true_pos / (true_pos + false_neg)

        stats_intermediate_df <- data.frame(precision = precision, recall = recall, method = method)

        stats_df <- rbind(stats_df, stats_intermediate_df)
    }

    stats_df$dataset = labels

    return(stats_df)
}

stats_df <- data.frame(precision = c(), recall = c(), method = c())

PAGA_df <- calculate_metrics(file_list, PAGA_prefix, PAGA_folder, method = "PAGA_umap") # nolint: line_length_linter.

monocle_df <- calculate_metrics(file_list, monocle_prefix, monocle_folder, method = "monocle_ddr") # nolint: line_length_linter.

slingshot_df <- calculate_metrics(file_list, Slingshot_prefix, Slingshot_folder, method = "slingshot_umap")

full_df <- rbind(stats_df, PAGA_df, monocle_df, slingshot_df)

full_df[is.na(full_df)] <- 0

plot_list = list()

#for each simulated dataset:
for (drop_rate in unique(full_df$dataset)){

    #isolate the precision and recall for just that dataset
    one_set <- subset(full_df, dataset == drop_rate)

    #make a plot of the precision and recall for each method
    plot <- ggplot(one_set, aes(x = precision, y = recall, color = method)) +
        geom_point(size = 5) +
        xlim(0, 0.5) +
        ylim(0, 0.5) +
        labs(x = "Precision", y = "Recall", title = "Precision vs Recall") +
        theme_minimal() +
        theme(
            plot.title = element_text(hjust = 0.5),
            legend.title = element_blank()
    )

    plot_list <- append(plot_list, list(plot))

}

#now we wrap the plots into one plot
all_scatter <- wrap_plots(plot_list, ncol = length(file_list)) & theme(legend.position = "right")

all_scatter <- all_scatter + plot_layout(guides = "collect")

#save monocle scatterplots
ggsave(paste0(out_dir, "/default_celltype_SvE3_zoom.png"), all_scatter, height = 5, width = 30)

#this works but I want to remove the legend and do it manually
