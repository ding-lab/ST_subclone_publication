#Basic commands to run Space ranger for spatial transcriptomics Data

#References: https://www.10xgenomics.com/support/software/space-ranger

#Reference Genome: https://www.10xgenomics.com/support/software/space-ranger/downloads#reference-downloads

#Sample_ID
sample=''

#Path to reference
reference='refdata-gex-GRCh38-2020-A'

# Probe set (FFPE): Human Probe Set v1
# Available at 10x Genomics website: https://www.10xgenomics.com/support/spatial-gene-expression-ffpe/documentation/workflows/cytassist-ffpe/steps/probe-sets/visium-ffpe-probe-sets-overview
Probe_set_file = "PATH/TO/Visium_Human_Transcriptome_Probe_Set_v1.0_GRCh38-2020-A.csv"

#Path to FASTQ files
fastq_dir='PATH/TO/FASTQ/'

# Path to output
output_folder="SPACERANGE_OUTPUT/FOLDER/${sample}"

# (Optional) in terms of sample with manual alignment
# Add 
# --loupe-alignment="PATH/TO/SAMPLE_manual_alignment.json" \

# OCT version
mkdir -p $output_folder && cd $output_folder && spaceranger count --id="$sample" \
--sample="$sample" \
--transcriptome="$reference" \
--fastqs="$fastq_dir" \
--image="PATH/TO/HE_IMAGE_FILE.tif" \
--slide="SLIDE_ID" \ # e.g. "V12U06-375"
--area="AREA_ID" \ # e.g. "B1"
--reorient-images \
--localcores=32 \
--localmem=150


# FFPE version
mkdir -p $output_folder && cd $output_folder && spaceranger count --id="$sample" \
--sample="$sample" \
--transcriptome="$reference" \
--probe-set="$Probe_set_file" \
--fastqs="$fastq_dir" \
--image="PATH/TO/HE_IMAGE_FILE.tif" \
--slide="SLIDE_ID" \ # e.g. "V12U06-375"
--area="AREA_ID" \ # e.g. "B1"
--reorient-images \
--localcores=32 \
--localmem=150

