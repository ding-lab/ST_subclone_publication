eval "$(conda shell.bash hook)"
conda activate seurat4.3
# Run all sample in 2_get_matching_spot with 1_get_multi_section_regions.r
PROJECT_FOLDER="/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/ST_subclone/30-PASTE2/1_PASTE2_result"
MATHING_OUT_PATH="$PROJECT_FOLDER/2_get_matching_spot/out"
ANALYSIS_PATH="$PROJECT_FOLDER/4_multi_section_regions_Ver20231031_N_min_3/"
SCRIPT_PATH="$ANALYSIS_PATH/script/1_get_multi_section_regions.r"

# Name format {piece_name}_min_dist.csv in MATHING_OUT_PATH folder
for file in $MATHING_OUT_PATH/*_min_dist.csv; do
    piece_name=$(basename "$file" _min_dist.csv)
    # do something with $piece_name
    echo $piece_name
    Rscript $SCRIPT_PATH \
    --analysis_folder $ANALYSIS_PATH \
    --piece_use $piece_name \
    --n_spot_minimum 3 \
    --matching_out_path $MATHING_OUT_PATH > $ANALYSIS_PATH/log/$piece_name.log 2> $ANALYSIS_PATH/log/$piece_name.err &
done
