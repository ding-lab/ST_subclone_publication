# Conda activate commot

## Load sample 
from optparse import OptionParser
parser = OptionParser()
parser.add_option("-f", "--file", dest="filename",
                  help="write report to FILE", metavar="FILE")
parser.add_option("-o", "--output", dest="output",
                  help="output folder", metavar="FILE", default = "./out/")
parser.add_option("-n", "--name", dest="name",
                  help="sample name", metavar="FILE", default = "untitled_sample")
parser.add_option("-s", "--species", dest="species",
                  help="speices to grab data from. human or mouse. Default human", metavar="FILE", default = "human")
parser.add_option("-t", "--dist_threshold", dest="dist_threshold",
                  help="Distance threshold used in spatial_communication function. Default = 1000", metavar="FILE", default = 1000)
parser.add_option("-d", "--database", dest="database",
                  help="Database to use to get LR paris. CellChat or CellPhoneDB_v4.0. Default CellChat", metavar="FILE", default = "CellChat")

                  


(options, args) = parser.parse_args()
print(options)

# From https://commot.readthedocs.io/en/latest/notebooks/visium-mouse_brain.html
import os
import gc
import ot
import pickle
import anndata
import scanpy as sc
import pandas as pd
import numpy as np
from scipy import sparse
from scipy.stats import spearmanr, pearsonr
from scipy.spatial import distance_matrix
import matplotlib.pyplot as plt
from pathlib import Path

import commot as ct
import squidpy as sq

# 0. parameter
root_out = f"{options.output}/{options.name}"
Path(root_out).mkdir(parents=True, exist_ok=True)
os.chdir(root_out)

# Load data
#adata = sc.datasets.visium_sge(sample_id='V1_Mouse_Brain_Sagittal_Posterior')
adata = sq.read.visium(f"{options.filename}")

# Preprocessing the data
adata.var_names_make_unique()
adata.raw = adata
sc.pp.normalize_total(adata, inplace=True)
sc.pp.log1p(adata)

adata_dis500 = adata.copy()

sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata = adata[:, adata.var.highly_variable]
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.4)
sc.pl.spatial(adata, color='leiden', save = "1_leiden.png")

# Spatial communication inference
# “Inference and analysis of cell-cell communication using CellChat.” Nature communications 12.1 (2021): 1-20.
df_cellchat = ct.pp.ligand_receptor_database(species=options.species, signaling_type=None, database=options.database) # set None to get all the LR pairs
print(df_cellchat.shape)

# also setup database name to use later
if options.database == "CellChat":
    database_name_use = 'cellchat'
elif options.database == "CellPhoneDB_v4.0":
    database_name_use = 'cellphonedb'

# Filter to keep only those express > 5%
df_cellchat_filtered = ct.pp.filter_lr_database(df_cellchat, adata_dis500, min_cell_pct=0.05)
print(df_cellchat_filtered.shape)
print(df_cellchat_filtered.head())

# spatial communiation inference. 
# First check if file exists
if os.path.isfile("./2_adata_LR_calculated.h5ad"):
    print("Loading spatial communication inference.. ")
    adata_dis500 = sc.read("./2_adata_LR_calculated.h5ad")
else:
    print("Starting spatial communication inference.. ")
    ct.tl.spatial_communication(adata_dis500,
        database_name=database_name_use, df_ligrec=df_cellchat_filtered, dis_thr=options.dist_threshold, heteromeric=True, pathway_sum=True) # change distance to 1000
    adata_dis500.write("./2_adata_LR_calculated.h5ad")

# Plot the spatial communication network
# Plot ALL The PATHWAYS
# pathways_to_plot = [a.replace("commot_cluster-leiden-cellchat-","") for a in adata_dis500.uns.keys() if 'commot_cluster-leiden-cellchat-' in a]
pathways_to_plot = set(df_cellchat_filtered.iloc[:,2].tolist())

for pathway_use in pathways_to_plot:
    print(pathway_use)
    # change dir
    Path(f"{root_out}/pathway/{pathway_use}").mkdir(parents=True, exist_ok=True)
    os.chdir(f"{root_out}/pathway/{pathway_use}")
    # first get the direction
    ct.tl.communication_direction(adata_dis500, database_name=database_name_use, pathway_name=pathway_use, k=5)
    # 3. then plot
    ct.pl.plot_cell_communication(adata_dis500, database_name=database_name_use, pathway_name=pathway_use, plot_method='grid', background_legend=True,
        scale=0.00003, ndsize=8, grid_density=0.4, summary='sender', background='image', clustering='leiden', cmap='Alphabet',
        normalize_v = True, normalize_v_quantile=0.995,
        filename = "3_sender.pdf")

    # Or summarize the signaling to the leiden clustering
    adata_dis500.obs['leiden'] = adata.obs['leiden']
    ct.tl.cluster_communication(adata_dis500, database_name=database_name_use, pathway_name=pathway_use, clustering='leiden',
        n_permutations=100)

    # 4. Plot with automatic node embedding
    ct.pl.plot_cluster_communication_network(adata_dis500, uns_names=[f'commot_cluster-leiden-cellchat-{pathway_use}'],
        nx_node_pos=None, nx_bg_pos=False, p_value_cutoff = 5e-2, filename=f'4_cluster.pdf', nx_node_cmap='Light24')
    # 5. Or with spatial node embedding.
    ct.tl.cluster_position(adata_dis500, clustering='leiden')
    ct.pl.plot_cluster_communication_network(adata_dis500, uns_names=[f'commot_cluster-leiden-cellchat-{pathway_use}'], clustering='leiden',
        nx_node_pos='cluster', nx_pos_idx=np.array([0, 1]), nx_bg_pos=True, nx_bg_ndsize=0.25, p_value_cutoff=5e-2,
        filename=f'5_cluster_spatial.pdf', nx_node_cmap='Light24')
# change dir back
os.chdir(root_out)


# # Downstream analysis
# # 6. Identify signaling DE genes
# # Van den Berge, Koen, et al. “Trajectory-based differential expression analysis for single-cell sequencing data.” Nature communications 11.1 (2020): 1-13.
# # AKA: TradeSeq
# adata_dis500 = sc.read_h5ad("./2_adata_LR_calculated.h5ad")
# adata = sq.read.visium(f"{options.filename}")
# adata_dis500.layers['counts'] = adata.X

# # Plot ALL available pathways
# available_pathways = set(df_cellchat_filtered.iloc[:,2].tolist())

# for pathway_use in available_pathways:
#     print(pathway_use)
#     # change dir
#     Path(f"{root_out}/{pathway_use}").mkdir(parents=True, exist_ok=True)
#     os.chdir(f"{root_out}/{pathway_use}")

#     # Look for genes that are differentially expressed with respect to PSAP signaling.
#     df_deg, df_yhat = ct.tl.communication_deg_detection(adata_dis500,
#         database_name = 'cellchat', pathway_name=pathway_use, summary = 'receiver')
#     import pickle
#     deg_result = {"df_deg": df_deg, "df_yhat": df_yhat}
#     with open(f'./6_deg_{pathway_use}.pkl', 'wb') as handle:
#         pickle.dump(deg_result, handle, protocol=pickle.HIGHEST_PROTOCOL)
#     # Cluster the downstream genes and visualize the expression trends with respect to increased level of received signal through PSAP pathway (horizontal axis).
#     with open(f"./6_deg_{pathway_use}.pkl", 'rb') as file:
#         deg_result = pickle.load(file)
#     df_deg_clus, df_yhat_clus = ct.tl.communication_deg_clustering(df_deg, df_yhat, deg_clustering_res=0.4)
#     top_de_genes_use = ct.pl.plot_communication_dependent_genes(df_deg_clus, df_yhat_clus, top_ngene_per_cluster=5,
#         filename=f'./6a_heatmap_deg_{pathway_use}.pdf', font_scale=1.2, return_genes=True)
#     # Plot some example signaling DE genes.
#     X_sc = adata_dis500.obsm['spatial']
#     fig, ax = plt.subplots(1,3, figsize=(15,4))
#     colors = adata_dis500.obsm['commot-cellchat-sum-receiver'][f'r-{pathway_use}'].values
#     idx = np.argsort(colors)
#     ax[0].scatter(X_sc[idx,0],X_sc[idx,1], c=colors[idx], cmap='coolwarm', s=10)
#     colors = adata_dis500[:,'Ctxn1'].X.toarray().flatten()
#     idx = np.argsort(colors)
#     ax[1].scatter(X_sc[idx,0],X_sc[idx,1], c=colors[idx], cmap='coolwarm', s=10)
#     colors = adata_dis500[:,'Gpr37'].X.toarray().flatten()
#     idx = np.argsort(colors)
#     ax[2].scatter(X_sc[idx,0],X_sc[idx,1], c=colors[idx], cmap='coolwarm', s=10)
#     ax[0].set_title('Amount of received signal')
#     ax[1].set_title('An example negative DE gene ()')
#     ax[2].set_title('An example positive DE gene ()')
#     # Further quantify impact of signaling on the DE genes
#     # The DE analysis above shows correlation between signaling and downstream gene expression. Now we further build random forest models with potential DE genes as targets and signaling as input and quantify the impact of signaling by feature importances.
#     df_impact_use = ct.tl.communication_impact(adata_dis500, database_name='cellchat', pathway_name = pathway_use,\
#         tree_combined = True, method = 'treebased_score', tree_ntrees=100, tree_repeat = 100, tree_method = 'rf', \
#         ds_genes = top_de_genes_use, bg_genes = 500, normalize=True)
#     # Visualize the impact scores as a heatmap. Here we only have two LR pairs in PSAP and we show the top 30 DE genes.
#     ct.pl.plot_communication_impact(df_impact_use, summary = 'receiver', top_ngene= 30, top_ncomm = 5, colormap='coolwarm',
#         font_scale=1.2, linewidth=0, show_gene_names=True, show_comm_names=True, cluster_knn=2,
#         filename = f'heatmap_impact_{pathway_use}.pdf')

# # change dir back
# os.chdir(root_out)

