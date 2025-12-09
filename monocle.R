#!/usr/bin/env Rscript
library(monocle)
library(igraph)
library(reshape)
library(profmem)

data_dir <- "/home/estagaman/benchmarking_project/data/trajectory/with_DE_genes/CDS"

#directory to save all files/results to
out_dir <- "/home/estagaman/benchmarking_project/test_simulated_data/monocle_ICA_mem_1end"

set.seed(123)

#check the first to make sure it's normal
#cds_object <- readRDS(paste0(data_dir, "/cds1.rds"))

#es.mef <- cds_object

#do quality control
#L <- log10(exprs(es.mef)+1)
#L[L==0] <- NA
#melted.dens.df <- reshape::melt(t(scale(t(L))))

#check_plot <- qplot(value, geom = 'density', data = melted.dens.df) + stat_function(fun = dnorm, size = 0.5, color = 'red') + xlab('Standardized log(Expression)') + ylab('Density')

#ggsave(paste0(out_dir, "/check_QC.png"), check_plot)

#if I need to do filtering, use this code:
#es.mef <- detectGenes(es.mef, min_expr = 10)

#expressed.genes <- rownames(fData(es.mef))[fData(es.mef)$num_cells_expressed >= 100]
#es.mef <- es.mef[expressed.genes, ]


#for running each individually

#i_options <- c("_glob1", "_glob2", "_glob3", "_glob4", "_glob5", "_glob6", "1", "2", "3", "4", "5", "6")

#es.mef <- readRDS(paste0(data_dir, "/cds", i_options[1], ".rds"))

#pseudotime_ordering(es.mef, i, out_dir)

#num = i_options[1]

#ok. now we create a loop for the pseudotime ordering
pseudotime_ordering <- function(es.mef, num, out_dir){

    fData(es.mef)$gene_short_name <- rownames(es.mef) #add column with gene name to feature data

    es.mef <- estimateSizeFactors(es.mef) #estimate library sizes

    es.mef <- reduceDimension(es.mef, reduction_method = "ICA") # Reduce dimensionality - this can take a long time

    #run order cells once
    ordered_1_start <- Sys.time()

    p_1 <- profmem({
        es.mef.test <- orderCells(es.mef, reverse = F)
    })

    ordered_1_end <- Sys.time()

    #to find the root state
    just_first_10 <- subset(pData(es.mef.test), pseudotime < 10)
    states_counted <- table(just_first_10$State)

    max_frequency_index <- which.max(states_counted)
    most_common_state <- names(states_counted[max_frequency_index])

    root_chosen = most_common_state

    ordered_2_start <- Sys.time()

    p_2 <- profmem({ #ordering cells while measuring memory allocation

        es.mef <- orderCells(es.mef.test, reverse = F, root_state = root_chosen, num_paths = 1)

    })

    ordered_2_end <- Sys.time()

    final_ordering_time <- (ordered_2_end - ordered_2_start) + (ordered_1_end - ordered_1_start)

    total_bytes <- max(c(sum(as.numeric(p_1$bytes), na.rm = TRUE), sum(as.numeric(p_2$bytes), na.rm = TRUE)))

    tree <- plot_cell_trajectory(es.mef, color_by = "pseudotime") # Plot trajectory, color it by the ground truth pseudotime

    ggsave(paste0(out_dir, "/tree_", num, ".png"), tree)

    #I want to extract the pseudotime assignments 
    write.csv(pData(es.mef), paste0(out_dir, "/results_", num,".csv"))

    #identify DE genes 
    DEG_results <- differentialGeneTest(es.mef, fullModelFormulaStr = "~Pseudotime",
        reducedModelFormulaStr = "~1", relative_expr = TRUE, cores = 1,
        verbose = FALSE)

    DE_str <- c()
    
    for (pval in DEG_results$pval){
        if (pval < 0.05){
            DE_str <- c(DE_str, "DE")
        } else {
            DE_str <- c(DE_str, "notDE")
        }
    }

    DEG_results$DE <- DE_str

    write.csv(DEG_results, paste0(out_dir, "/DEG_", num,".csv"))

    #at the end, save all the time and space complexity info for each run
    stats_df <- data.frame("ordering_time" = final_ordering_time, "memory" = total_bytes)

    write.csv(stats_df, paste0(out_dir, "/stats_", num, ".csv"))
}

#run pseudotime ordering for each simulated file, saving results to out_dir
for (i in c("glob1", "glob2", "glob3", "glob4", "glob5", "glob6", "1", "2", "3", "4", "5", "6")){
    es.mef <- readRDS(paste0(data_dir, "/cds_", i, ".rds"))

    pseudotime_ordering(es.mef, i, out_dir)
}
