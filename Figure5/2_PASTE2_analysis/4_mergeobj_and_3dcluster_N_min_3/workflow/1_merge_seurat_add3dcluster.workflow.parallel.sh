eval "$(conda shell.bash hook)"
conda activate seurat4.3

# Name format {piece_name}_min_dist.csv in MATHING_OUT_PATH folder
export PROJECT_FOLDER=""
export MATHING_OUT_PATH="$PROJECT_FOLDER/2_get_matching_spot/out"
export VOLUMN_CLUSTER_PATH="$PROJECT_FOLDER/3_multi_section_regions_Ver20231031_N_min_3/out"
export MERGED_OBJECT_OUT="$PROJECT_FOLDER/4_mergeobj_and_3dcluster_N_min_3/"
export SCRIPT_PATH="$MERGED_OBJECT_OUT/script/1_merge_seurat_add3dcluster.r"

ls $MATHING_OUT_PATH/*_min_dist.csv | parallel --jobs 4 \
    'piece_name=$(basename {} _min_dist.csv); \
    echo $piece_name; \
    Rscript $SCRIPT_PATH \
    --analysis_folder $MERGED_OBJECT_OUT \
    --piece_use $piece_name \
    --matching_out_path $MATHING_OUT_PATH \
    --volumncluster_path $VOLUMN_CLUSTER_PATH > $MERGED_OBJECT_OUT/log/$piece_name.log 2> $MERGED_OBJECT_OUT/log/$piece_name.err'

# 2. Make plot
PLOT_SCRIPT_PATH="$PROJECT_FOLDER/5_mergeobj_and_3dcluster/script/2_plot_3d_cluster.r"
export PLOT_SCRIPT_PATH
export MERGED_OBJECT_OUT
export MATHING_OUT_PATH

ls $MATHING_OUT_PATH/*_min_dist.csv | parallel -j 4 '
    piece_name=$(basename {} _min_dist.csv)
    echo $piece_name
    Rscript $PLOT_SCRIPT_PATH \
    --analysis_folder $MERGED_OBJECT_OUT \
    --piece_use $piece_name \
    --merged_obj_path $MERGED_OBJECT_OUT >> $MERGED_OBJECT_OUT/log/$piece_name.log 2>> $MERGED_OBJECT_OUT/log/$piece_name.err
'


