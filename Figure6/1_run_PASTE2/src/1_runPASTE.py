import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import scanpy as sc
import paste as pst
import pandas as pd
import squidpy as sq


# use optparse to parse arguments
import optparse
parser = optparse.OptionParser()
# options for slide1 and slide2 paths
parser.add_option('-1', '--slide1_path', action="store", dest="slide1_path", help="slide1 path")
parser.add_option('-2', '--slide2_path', action="store", dest="slide2_path", help="slide2 path")
# options for slide1 and slide2 names
parser.add_option('-a', '--slide1_name', action="store", dest="slide1_name", help="slide1 name")
parser.add_option('-b', '--slide2_name', action="store", dest="slide2_name", help="slide2 name")
# options for opt path
parser.add_option('-o', '--out_root', action="store", dest="out_root", help="output path")
# option for paste2 path, default to "/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/ST_subclone/30-PASTE2/3_PASTE2_runs/0_paste2_github/src/"
parser.add_option('-p', '--paste2_path', action="store", 
                    default="/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/ST_subclone/30-PASTE2/3_PASTE2_runs/0_paste2_github/src/",
                    dest="paste2_path", help="paste2 path")
# Option to select which part to run. Default only PASTE1. 
parser.add_option('-r', '--run_part', action="store", default="123", dest="run_part", help="run part. 1=PASTE1, 2=PASTE2, 3=PASTE2-histology. 12=PASTE1+PASTE2. 13=PASTE1+PASTE2-histology. 23=PASTE2+PASTE2-histology. 123=PASTE1+PASTE2+PASTE2-histology")
# Set maxinum number of iterations for PASTE1
parser.add_option('-m', '--max_iter', action="store", default=1000, dest="max_iter", help="max number of iterations for PASTE1")

options, args = parser.parse_args()
slide1_path = options.slide1_path
slide2_path = options.slide2_path
slide1_name = options.slide1_name
slide2_name = options.slide2_name
out_root = options.out_root
paste2_path = options.paste2_path
run_part = options.run_part
max_iter = int(options.max_iter)

# create sample name
sample_name_list = [slide1_name, slide2_name]
sample_name_list.sort()
sample_run_name = f"{sample_name_list[0]}_{sample_name_list[1]}"
print(f"Running PASTE for {sample_run_name}")

# create out_path
out_path = f"{out_root}/{sample_run_name}"
import os
if not os.path.exists(out_path):
    os.makedirs(out_path)

# Load and Preprocess slices
def load_and_preprocess(sample_path: str):
    visium_sample = sq.read.visium(sample_path)
    sc.pp.filter_genes(visium_sample, min_counts = 15)
    sc.pp.filter_cells(visium_sample, min_counts = 100)
    # make variable name unique
    visium_sample.var_names_make_unique()
    return visium_sample


# Load visium
file_list = [
    slide1_path,
    slide2_path
]

slices = [load_and_preprocess(file_path) for file_path in file_list]
slice1 = slices[0]
slice2 = slices[1]

# Arguments from command line PASTE: https://github.com/raphael-group/paste/blob/main/paste-cmd-line.py
args_paste={
    'alpha' :0.1,
    'cost' :"kl",
}

# Plot parameters
spot_size = 70

##############################################
# Use PASTE2
##############################################
# preprocess
# first add rgb to slices
slice1.obsm['rgb'] = slice1.obsm['spatial']
slice2.obsm['rgb'] = slice2.obsm['spatial']

# import paste2
# add paste2 to the path
import sys
sys.path.append(paste2_path)
from paste2.PASTE2 import partial_pairwise_align, partial_pairwise_align_histology

"""
    Optimal partial alignment of two slices using both gene expression and histological image information.

    sliceA, sliceB must be AnnData objects that contain .obsm['rgb'], which stores the RGB value of each spot in the histology image.
"""
if(run_part.find("2") != -1):
    # Check if output file exists already
    output_part2_path = f"{out_path}/2_PASTE2_partial_pairwise_numpy.npy"
    if os.path.exists(output_part2_path):
        print("PASTE2 already run, loading results")
        slice12_al = np.load(output_part2_path)
    else:
        # Run PASTE2 alignment
        print("Running pairwise alignment from PASTE2")
        slice12_al = partial_pairwise_align(slice1, slice2, s= args_paste['alpha'], dissimilarity=args_paste['cost']) # this will take a while
        # save both results as npy
        np.save(output_part2_path, slice12_al)

    # PLOT ------------------------------------------------
    # To visualize the alignment you can stack the slices 
    # according to the alignment pi
    slices, pis = [slices[0], slices[1]], [slice12_al]
    new_slices = pst.stack_slices_pairwise(slices, pis)

    slice_colors = ['#e41a1c','#377eb8']
    plt.figure(figsize=(7,7))
    for i in range(len(new_slices)):
        pst.plot_slice(new_slices[i],slice_colors[i],s=spot_size)

    plt.legend(handles=[mpatches.Patch(color=slice_colors[0], label=sample_name_list[0]),
                        mpatches.Patch(color=slice_colors[1], label=sample_name_list[1])])
    plt.gca().invert_yaxis()
    plt.axis('off')
    #plt.show()
    plt.savefig(f"{out_path}/2_PASTE2_partial_pairwise.pdf", format="pdf", bbox_inches="tight")

    # Generate new coordinates
    # From the command line script: https://github.com/raphael-group/paste/blob/main/paste-cmd-line.py
    n_slices = 2
    for i in range(n_slices):
        # Format to match original tissue_positions.csv
        # barcode,in_tissue,array_row,array_col,pxl_row_in_fullres,pxl_col_in_fullres,pxl_row_in_fullres_new,pxl_col_in_fullres_new
        df_tissue = slices[i].obs[['in_tissue','array_row','array_col']] #.reset_index().rename(columns = {'index':'barcode'})
        df_coord_orig = pd.DataFrame(slices[i].obsm['spatial'], index = df_tissue.index, columns = ['pxl_row_in_fullres','pxl_col_in_fullres'])
        df_coord_new  = pd.DataFrame(new_slices[i].obsm['spatial'], index = df_tissue.index, columns = ['pxl_row_in_fullres_new','pxl_col_in_fullres_new'])
        df_coord = pd.concat([df_tissue,df_coord_orig, df_coord_new], axis = 1).reset_index().rename(columns = {'index':'barcode'})
        df_coord.to_csv(f"{out_path}/2_PASTE2_partial_pairwise_" + sample_name_list[i] + "_new_tissue_positions.csv", index = False)



##############################################
# PASTE2 algnment using histology
##############################################
if(run_part.find("3") != -1):
    # Check if output file exists already
    npy_out_part3_path=f"{out_path}/3_PASTE2_partial_pairwise_histology_numpy.npy"
    if(os.path.isfile(npy_out_part3_path)):
        print("PASTE2 histology already run, loading results")
        # Load file directly# Load PASTE2 histology alignment
        slice12_al_hist = np.load(npy_out_part3_path, allow_pickle=True)
    else:
        # Run PASTE2 histology alignment
        print("Running pairwise alignment from PASTE2 using histology")
        slice12_al_hist = partial_pairwise_align_histology(slice1, slice2, s= args_paste['alpha'], dissimilarity=args_paste['cost']) # this will take a while
        np.save(npy_out_part3_path, slice12_al_hist)

    
    # Plot ------------------------------------------------
    # To visualize the alignment you can stack the slices
    # according to the alignment pi
    slices, pis = [slices[0], slices[1]], [slice12_al_hist]
    new_slices = pst.stack_slices_pairwise(slices, pis)

    slice_colors = ['#e41a1c','#377eb8']
    plt.figure(figsize=(7,7))
    for i in range(len(new_slices)):
        pst.plot_slice(new_slices[i],slice_colors[i],s=spot_size)

    plt.legend(handles=[mpatches.Patch(color=slice_colors[0], label=sample_name_list[0]),
                        mpatches.Patch(color=slice_colors[1], label=sample_name_list[1])])
    plt.gca().invert_yaxis()
    plt.axis('off')
    plt.savefig(f"{out_path}/3_PASTE2_partial_pairwise_histology.pdf", format="pdf", bbox_inches="tight")

    n_slices = 2
    for i in range(n_slices):
        # Format to match original tissue_positions.csv
        # barcode,in_tissue,array_row,array_col,pxl_row_in_fullres,pxl_col_in_fullres,pxl_row_in_fullres_new,pxl_col_in_fullres_new
        df_tissue = slices[i].obs[['in_tissue','array_row','array_col']] #.reset_index().rename(columns = {'index':'barcode'})
        df_coord_orig = pd.DataFrame(slices[i].obsm['spatial'], index = df_tissue.index, columns = ['pxl_row_in_fullres','pxl_col_in_fullres'])
        df_coord_new  = pd.DataFrame(new_slices[i].obsm['spatial'], index = df_tissue.index, columns = ['pxl_row_in_fullres_new','pxl_col_in_fullres_new'])
        df_coord = pd.concat([df_tissue,df_coord_orig, df_coord_new], axis = 1).reset_index().rename(columns = {'index':'barcode'})
        df_coord.to_csv(f"{out_path}/3_PASTE2_partial_pairwise_histology_" + sample_name_list[i] + "_new_tissue_positions.csv", index = False)
