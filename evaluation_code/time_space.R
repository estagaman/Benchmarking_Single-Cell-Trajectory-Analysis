#plotting time and space complexity
library(patchwork)
library(ggplot2)

#folder with monocle time and space statistics, with beginning prefix of all filenames with statistics
monocle_folder <- "/home/estagaman/benchmarking_project/test_simulated_data/monocle_DDR"
monocle_prefix = "stats_"

#folder with PAGA time and space statistics, with beginning prefix of all filenames with statistics
PAGA_folder <- "/home/estagaman/benchmarking_project/test_simulated_data/PAGA_out_umap_30/seurat"
PAGA_prefix = "stats_"

#folder with slingshot time and space statistics, with beginning prefix of all filenames with statistics
Slingshot_folder <- "/home/estagaman/benchmarking_project/test_simulated_data/slingshot_out_umap"
Slingshot_prefix = "stats_"

#do we want to look at global dropouts or cell-type-specific dropouts?
check <- "global"

#specify output folder
out_dir <- "/home/estagaman/benchmarking_project/test_simulated_data/time_space_results"

#picks filenames and dataset settings labels based off whether we checked global or cell_type
if (check == "global"){

    file_list <- c("1", "glob1", "glob2", "glob3", "glob4", "glob5", "glob6")
    labels <- c("no dropout", "0.25 dropout", "0.5 dropout", "0.75 dropout", "1 dropout", "1.25 dropout", "1.5 dropout")

} else if (check == "cell_type"){

    file_list <- c("1", "2", "3", "4", "5", "6")
    labels <- c("no dropout", "[0.1, 0.15, 0.2, 0.25]", "[0.5, 1, 1.25, 1.5]", "[1, 1.25, 1.5, 1.75]", "[1.5, 1.75, 2, 2.25]", "[0.1, 0.3, 0.5, 0.7]")

}

#calculate time and space metrics
calculate_metrics <- function(file_list, prefix, results_folder, method){

    #initialize dataset for statistics
    stats_df <- data.frame(time = c(), memory = c(), method = c())

    #for each file:
    for (file in file_list){

        #read in the results file
        res <- read.csv(paste0(results_folder, "/", prefix, file, ".csv"))

        #check for ordering time
        ordering_time <- res$ordering_time[1] #only one time should be reported

        #check for ordering memory
        ordering_memory <- res$memory[1]

        #collect these stats in a data frame
        stats_intermediate_df <- data.frame(time = ordering_time, memory = ordering_memory, method = method)

        #bind them together with stats from other files
        stats_df <- rbind(stats_df, stats_intermediate_df)
    }

    #add dataset dropout labels for description
    stats_df$dataset = labels

    #return the statistics
    return(stats_df)
}

#calculate PAGA metrics
PAGA_df <- calculate_metrics(file_list, PAGA_prefix, PAGA_folder, method = "PAGA_UMAP") # nolint: line_length_linter.

#calculate Monocle metrics
monocle_df <- calculate_metrics(file_list, monocle_prefix, monocle_folder, method = "monocle_DDR") # nolint: line_length_linter.

#calculate Slingshot metrics
slingshot_df <- calculate_metrics(file_list, Slingshot_prefix, Slingshot_folder, method = "slingshot_UMAP")

#combine dfs together 
full_df <- rbind(PAGA_df, monocle_df, slingshot_df)

#make plot of the memory: 
memory_plot <- ggplot(full_df, aes(x = factor(dataset, levels = labels), y = memory, color = method)) +
  geom_line(size = 1) + #line between points
  geom_point(size = 2) + #plot points
  labs(title = "Memory Usage vs Dropout", #title
    x = "Dropout rate", #x and y axis labels
    y = "Memory (bytes)") +
  theme_classic()

#make plot of time
time_plot <- ggplot(full_df, aes(x = factor(dataset, levels = labels), y = time, color = method)) +
  geom_line(size = 1) + #line between points
  geom_point(size = 2) + #plot points
  labs(title = "Time vs Dropout", #title
    x = "Dropout rate", #x and y axis labels
    y = "Time (seconds)") +
  theme_classic()

#save the plot
ggsave(paste0(out_dir, "/time_space_", check, "_default.png"), memory_plot + time_plot, width = 20, height = 5)

#save the statistics in a csv
write.csv(full_df, paste0(out_dir, "/time_space_", check, "_default.csv")
