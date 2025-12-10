#!/usr/bin/env Rscript

#load necessary packages
library(monocle)
library(igraph)
library(reshape)
library(profmem)

#directory where CellDataSets are stored with cells you would like to order
data_dir <- "/home/estagaman/benchmarking_project/data/trajectory/with_DE_genes/CDS"

#directory to save all files/results to
out_dir <- "/home/estagaman/benchmarking_project/test_simulated_data/monocle_DDR"

#set a seed for reproducibility
set.seed(123)

#ok. now we create a loop for the pseudotime ordering
pseudotime_ordering <- function(es.mef, num, out_dir){ #takes in cell data set, the identifier of the specific file we're analyzing, and output directory

    #add column with gene name to feature data
    fData(es.mef)$gene_short_name <- rownames(es.mef) 

    #estimate library sizes
    es.mef <- estimateSizeFactors(es.mef) #estimate library sizes

    # Reduce dimensionality - this can take a long time
    es.mef <- reduceDimension(es.mef, reduction_method = "DDRTree") 

    #run order cells once
    ordered_1_start <- Sys.time() #start timer

    p_1 <- profmem({ #use profmem to record memory usage
        es.mef.test <- orderCells(es.mef, reverse = F)
    })

    ordered_1_end <- Sys.time() #end timer

    #to find the root state
    just_first_10 <- subset(pData(es.mef.test), pseudotime < 10) #isolate the cells from the first 10% of true pseudotime
    states_counted <- table(just_first_10$State) #find the clusters/states that contain those cells

    max_frequency_index <- which.max(states_counted) #find the cluster containing the highest frequency of the selected cells
    most_common_state <- names(states_counted[max_frequency_index]) #this is our root cluster

    root_chosen = most_common_state #this is our root cluster

    ordered_2_start <- Sys.time() #start timer again

    p_2 <- profmem({ #ordering cells while measuring memory allocation

        #order cells again, but specify the root state to begin trajectory from
        es.mef <- orderCells(es.mef.test, reverse = F, root_state = root_chosen)

    })

    ordered_2_end <- Sys.time() #stop the timer

    #compute total ordering time
    final_ordering_time <- (ordered_2_end - ordered_2_start) + (ordered_1_end - ordered_1_start)

    #compute maximum memory usage between the two runs of orderCells()
    total_bytes <- max(c(sum(as.numeric(p_1$bytes), na.rm = TRUE), sum(as.numeric(p_2$bytes), na.rm = TRUE)))

    #plot the cell trajectory
    tree <- plot_cell_trajectory(es.mef, color_by = "pseudotime") # Plot trajectory, color it by the ground truth pseudotime

    #save the plot if you want to view the trajectory
    ggsave(paste0(out_dir, "/tree_", num, ".png"), tree)

    #I want to extract the pseudotime assignments and save them to file pseudotime_ + identifier of the cell data set
    write.csv(pData(es.mef), paste0(out_dir, "/pseudotime_", num,".csv"))

    #at the end, save all the time and space complexity info for each run
    stats_df <- data.frame("ordering_time" = final_ordering_time, "memory" = total_bytes)

    #write to a csv
    write.csv(stats_df, paste0(out_dir, "/stats_", num, ".csv"))
}

#run pseudotime ordering for each simulated file, saving results to out_dir
#this vector contains the unique part of each filename we want to test. If your files are named differentially than my simulated data, you would need to update these
for (i in c("glob1", "glob2", "glob3", "glob4", "glob5", "glob6", "1", "2", "3", "4", "5", "6")){

    #we know all files start in cds_ and end in .rds, so we just need to specify the middle section
    es.mef <- readRDS(paste0(data_dir, "/cds_", i, ".rds")) #for example, if i is glob5, the filename accessed would be cds_glob5.rds. Change this to match your filename format

    #do the pseudotime ordering using function above
    pseudotime_ordering(es.mef, i, out_dir)
}
