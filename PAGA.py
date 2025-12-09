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

out_dir = "/home/estagaman/benchmarking_project/test_simulated_data/PAGA_out_ica_30"

file_list = ["glob1", "glob2", "glob3", "glob4", "glob5", "glob6", "1", "2", "3", "4", "5", "6"]

for i in file_list:
    
    #read in the data
    adata = sc.read_h5ad("/home/estagaman/benchmarking_project/data/trajectory/with_DE_genes/h5ad/seurat_" + i + ".h5ad")

    #identify the seurat cluster holding the beginning of pseudotime
    just_first_10 = adata.obs['pseudotime'] < 5
    filt_first_10 = adata.obs[just_first_10]
    most_common_cluster = filt_first_10['seurat_clusters'].mode().iloc[0]
    root_chosen = str(most_common_cluster)

    # Compute neighbors using ICA
    sc.pp.neighbors(adata, use_rep='X_ica', n_neighbors=30)  # you can adjust n_neighbors

    # Now you can run draw_graph
    adata.obs['seurat_clusters'] = adata.obs['seurat_clusters'].astype('category') #converting seurat clusters to categorical
    sc.tl.draw_graph(adata)  # no init_pos
    sc.pl.draw_graph(adata, color=['seurat_clusters','group', 'pseudotime'], legend_loc='on data', legend_fontsize='xx-small')

    #perform Louvain clustering
    sc.tl.leiden(adata, resolution=1.0)

    #convert leiden clusters to categorical and identify root cluster
    adata.obs['leiden'] = adata.obs['leiden'].astype(str).astype('category')
    just_first_10 = adata.obs['pseudotime'] < 5
    filt_first_10 = adata.obs[just_first_10]
    most_common_leiden = filt_first_10['leiden'].mode().iloc[0]
    root_chosen_leiden = str(most_common_leiden)

    #run PAGA with leiden results - this is where I want to measure time and space complexity

    leiden_start_time = time.perf_counter()

    tracemalloc.start() #measuring space complexity

    sc.tl.paga(adata, groups='leiden')
    #sc.pl.paga(adata, labels="leiden")

    #give it a root state
    adata.uns['iroot'] = np.flatnonzero(adata.obs['leiden']  == root_chosen_leiden)[0]

    #do diffusion pseudotime for exact pseudotime value
    sc.tl.dpt(adata)

    snapshot = tracemalloc.take_snapshot()
    tracemalloc.stop() #stop measuring space complexity

    leiden_end_time = time.perf_counter()

    #calculate memory complexity
    stats = snapshot.statistics('lineno')

    total_memory_leiden = sum([stat.size for stat in stats])

    #visualize this
    sc.pl.draw_graph(adata, color=['seurat_clusters', 'pseudotime', 'dpt_pseudotime'], legend_loc='on data', legend_fontsize= 'x-small')

    #save the pseudotime values assigned
    pseudotime_values = adata.obs[["pseudotime", "dpt_pseudotime"]]

    leiden_file = os.path.join(out_dir, "leiden/pseudotime" + i + ".csv")
    pseudotime_values.to_csv(leiden_file, index=True)

    #####
    #run PAGA again but with seurat clusters instead

    seurat_start_time = time.perf_counter()

    tracemalloc.start() #measuring space complexity

    sc.tl.paga(adata, groups='seurat_clusters')

    #give it a root state
    adata.uns['iroot'] = np.flatnonzero(adata.obs['seurat_clusters']  == int(root_chosen))[0]

    #do diffusion pseudotime for exact pseudotime value
    sc.tl.dpt(adata)

    snapshot = tracemalloc.take_snapshot()
    tracemalloc.stop() #stop measuring space complexity

    seurat_end_time = time.perf_counter()

    #calculate memory used
    stats = snapshot.statistics('lineno')
    total_memory_seurat = sum([stat.size for stat in stats])

    #visualize this
    sc.pl.draw_graph(adata, color=['seurat_clusters', 'pseudotime', 'dpt_pseudotime'], legend_loc='on data', legend_fontsize= 'x-small')

    #save the pseudotime values assigned
    pseudotime_values = adata.obs[["pseudotime", "dpt_pseudotime"]]

    seurat_file = os.path.join(out_dir, "seurat/pseudotime" + i + ".csv")
    pseudotime_values.to_csv(seurat_file, index=True)

    #save time and space complexity
    leiden_time = leiden_end_time - leiden_start_time
    seurat_time = seurat_end_time - seurat_start_time 

    time_space_leiden = pd.DataFrame({"ordering_time": [leiden_time], "memory":[total_memory_leiden]})
    time_space_seurat = pd.DataFrame({"ordering_time": [seurat_time], "memory": [total_memory_seurat]})

    leiden_path = os.path.join(out_dir, "leiden/stats_" + i + ".csv")
    seurat_path = os.path.join(out_dir, "seurat/stats_" + i + ".csv")

    time_space_leiden.to_csv(leiden_path)
    time_space_seurat.to_csv(seurat_path)
