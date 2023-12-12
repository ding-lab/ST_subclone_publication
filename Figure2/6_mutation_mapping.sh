#!/bin/bash
# 10Xmapping run for Visium with matching WES somatic mutation calls
# Jingxian Clara Liu
# https://github.com/ding-lab/10Xmapping

# Start interactive session (Compute1)
# bsub -g compute-dinglab -q dinglab-interactive -is -a 'docker(scao/dailybox)' /bin/bash
# bsub -G compute-dinglab -q general-interactive -Is -a 'docker(cliu97/10xmapping)' /bin/bash

# env: genomics (katmai)
eval "$(conda shell.bash hook)"
conda activate genomics

script_dir="/diskmnt/Projects/Users/cliu/software/10Xmapping/"
work_dir=$1

# column 16: Tumor_Sample_Barcode
# maf=${work_dir}/03.20221129_ST.dnp.annotated.maf
# slice_id="HT260C1"
# sample_id="HT260C1"
# bam="/diskmnt/Datasets/Spatial_Transcriptomics/outputs_OCT/Human/HT260C1/H1/HT260C1-Th1K1Fc2U1Z1Bs1/outs/possorted_genome_bam.bam"
maf=$2
slice_id=$3
sample_id=$4
bam=$5

mkdir -p ${work_dir}/${slice_id}

cd ${script_dir}
# step 1: extracting reads contains the reference allele and variant allele of a somatic variant
tmp_maf=${work_dir}/${slice_id}/${slice_id}_filtered.maf
# head -n 1 ${maf} > ${tmp_maf}
# awk -v sample="${sample_id}" '$16 ~ sample { print }' ${maf} >> ${tmp_maf}
# perl 10Xmapping.pl --bam ${bam} --maf ${tmp_maf} --out ${work_dir}/${slice_id}/${slice_id}_step1.out --mapq 255
perl 10Xmapping_clara.pl --bam ${bam} --maf ${tmp_maf} --out ${work_dir}/${slice_id}/${slice_id}_step1.out --mapq 255

# step 2: get the corresponding table between ref allele, var allele and 10X barcode
perl parse_scrna_bc.pl ${work_dir}/${slice_id}/${slice_id}_step1.out ${work_dir}/${slice_id}/${slice_id}_step2.out

# step 3: output number of reads supporting reference and variant alleles and VAF
# perl rc.pl ${work_dir}/${slice_id}/${slice_id}_step2.out ${work_dir}/${slice_id}/${slice_id}_step3.out
perl rc_clara.pl ${work_dir}/${slice_id}/${slice_id}_step2.out ${work_dir}/${slice_id}/${slice_id}_step3.out

# step 4: barcode -> variant mapping
perl format_rc_clara.pl ${work_dir}/${slice_id}/${slice_id}_step2.out ${work_dir}/${slice_id}/${slice_id}_step4.out
