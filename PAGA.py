#import all necessary libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as pl
from matplotlib import rcParams
import scanpy as sc
import scipy
import numpy as np
import matplotlib.pyplot as plt
import warnings
import time
import os
import tracemalloc

#specify your output directory name
out_dir = "/home/estagaman/benchmarking_project/test_simulated_data/PAGA_out_ica_30"

#specify the unique element in each of your filenames
#for example, my simulated data filenames are formatted like so: "seurat_glob1.h5ad", "seurat_glob2.h5ad", "seurat_5.h5ad" 
#in this case I only include the section between the underscore and dot, which is unique to each file
file_list = ["glob1", "glob2", "glob3", "glob4", "glob5", "glob6", "1", "2", "3", "4", "5", "6"]

#for each file: 
for i in file_list:
    
    #read in the data
    adata = sc.read_h5ad("/home/estagaman/benchmarking_project/data/trajectory/with_DE_genes/h5ad/seurat_" + i + ".h5ad")

    #identify the seurat cluster holding the beginning of pseudotime
    just_first_10 = adata.obs['pseudotime'] < 5 #look at cells with true pseudotime less than 5
    filt_first_10 = adata.obs[just_first_10] #find the rows matching these cells
    most_common_cluster = filt_first_10['seurat_clusters'].mode().iloc[0] #find the most common seurat cluster these cells belong to
    root_chosen = str(most_common_cluster) #use this as the root cluster

    # Compute neighbors using UMAP
    sc.pp.neighbors(adata, use_rep='X_umap', n_neighbors=30)  # you can adjust n_neighbors if you want

    # Now you can run draw_graph
    adata.obs['seurat_clusters'] = adata.obs['seurat_clusters'].astype('category') #converting seurat clusters to categorical
    sc.tl.draw_graph(adata)  # draw the graph without any initial root state
    sc.pl.draw_graph(adata, color=['seurat_clusters','group', 'pseudotime'], legend_loc='on data', legend_fontsize='xx-small') #draw the graph and color code by seurat cluster, cell group, and true pseudotime for inspection

    #start timer
    seurat_start_time = time.perf_counter()

    #start measuring space complexity
    tracemalloc.start()

    #draw the PAGA graph using the seurat clusters as groups
    sc.tl.paga(adata, groups='seurat_clusters')

    #give it our chosen root state
    adata.uns['iroot'] = np.flatnonzero(adata.obs['seurat_clusters']  == int(root_chosen))[0]

    #do diffusion pseudotime for exact pseudotime value
    sc.tl.dpt(adata)

    #get the memory complexity
    snapshot = tracemalloc.take_snapshot()
    tracemalloc.stop() #stop measuring memory complexity

    #stop the timer
    seurat_end_time = time.perf_counter()

    #calculate memory used
    stats = snapshot.statistics('lineno')
    total_memory_seurat = sum([stat.size for stat in stats])

    #visualize these pseudotime results so I can check if they make sense
    sc.pl.draw_graph(adata, color=['seurat_clusters', 'pseudotime', 'dpt_pseudotime'], legend_loc='on data', legend_fontsize= 'x-small')

    #save the pseudotime values assigned
    pseudotime_values = adata.obs[["pseudotime", "dpt_pseudotime"]]

    #save the assigned pseudotime values to a csv with column "dpt_pseudotime"
    seurat_file = os.path.join(out_dir, "pseudotime_" + i + ".csv")
    pseudotime_values.to_csv(seurat_file, index=True)

    #save time and space complexity
    seurat_time = seurat_end_time - seurat_start_time 
    time_space_seurat = pd.DataFrame({"ordering_time": [seurat_time], "memory": [total_memory_seurat]})
    seurat_path = os.path.join(out_dir, "stats_" + i + ".csv")
    time_space_seurat.to_csv(seurat_path)
