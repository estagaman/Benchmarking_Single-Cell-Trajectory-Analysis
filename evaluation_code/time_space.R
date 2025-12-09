#plotting time and space complexity
library(patchwork)
library(ggplot2)

#Monocle
monocle_folder <- "/home/estagaman/benchmarking_project/test_simulated_data/monocle_DDR_mem"
monocle_prefix = "stats_"

#PAGA
PAGA_folder <- "/home/estagaman/benchmarking_project/test_simulated_data/PAGA_out_umap_30/seurat"
PAGA_prefix = "stats_"

#Slingshot
Slingshot_folder <- "/home/estagaman/benchmarking_project/test_simulated_data/slingshot_out_umap"
Slingshot_prefix = "stats_"

#do we want to look at global dropouts or cell-type-specific dropouts?
check <- "global"

#specify output folder
out_dir <- "/home/estagaman/benchmarking_project/test_simulated_data/time_space_results"

if (check == "global"){

    file_list <- c("1", "glob1", "glob2", "glob3", "glob4", "glob5", "glob6")
    labels <- c("no dropout", "0.25 dropout", "0.5 dropout", "0.75 dropout", "1 dropout", "1.25 dropout", "1.5 dropout")

} else if (check == "cell_type"){

    file_list <- c("1", "2", "3", "4", "5", "6")
    labels <- c("no dropout", "[0.1, 0.15, 0.2, 0.25]", "[0.5, 1, 1.25, 1.5]", "[1, 1.25, 1.5, 1.75]", "[1.5, 1.75, 2, 2.25]", "[0.1, 0.3, 0.5, 0.7]")

}

calculate_metrics <- function(file_list, prefix, results_folder, method){

    stats_df <- data.frame(time = c(), memory = c(), method = c())

    for (file in file_list){

        #read in the results file
        res <- read.csv(paste0(results_folder, "/", prefix, file, ".csv"))

        ordering_time <- res$ordering_time[1] #only one time should be reported

        ordering_memory <- res$memory[1]

        stats_intermediate_df <- data.frame(time = ordering_time, memory = ordering_memory, method = method)

        stats_df <- rbind(stats_df, stats_intermediate_df)
    }

    stats_df$dataset = labels

    return(stats_df)
}

PAGA_df <- calculate_metrics(file_list, PAGA_prefix, PAGA_folder, method = "PAGA_ICA") # nolint: line_length_linter.

monocle_df <- calculate_metrics(file_list, monocle_prefix, monocle_folder, method = "monocle_ICA") # nolint: line_length_linter.

slingshot_df <- calculate_metrics(file_list, Slingshot_prefix, Slingshot_folder, method = "slingshot_ICA")

full_df <- rbind(PAGA_df, monocle_df, slingshot_df)

#make plot of the memory: 
memory_plot <- ggplot(full_df, aes(x = factor(dataset, levels = labels), y = memory, color = method)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(title = "Memory Usage vs Dropout",
    x = "Dropout rate",
    y = "Memory (bytes)") +
  theme_classic()

time_plot <- ggplot(full_df, aes(x = factor(dataset, levels = labels), y = time, color = method)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(title = "Time vs Dropout",
    x = "Dropout rate",
    y = "Time (seconds)") +
  theme_classic()

ggsave(paste0(out_dir, "/celltype_default.png"), memory_plot + time_plot, width = 20, height = 5)
