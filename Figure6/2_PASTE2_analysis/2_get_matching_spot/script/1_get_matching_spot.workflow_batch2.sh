eval "(conda shell.bash hook)" # activate conda environment
conda activate seurat4.3

PROJECT_PATH=""
analysis_path="$PROJECT_PATH/2_get_matching_spot/script"

mkdir -p $analysis_path/log
Rscript --no-save --no-restore --verbose $analysis_path/1_get_matching_spot_batch2.r > $analysis_path/log/run.log 2> $analysis_path/log/run.err