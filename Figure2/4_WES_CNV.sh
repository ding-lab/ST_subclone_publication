#!/bin/bash
# Run WES processing pipeline including alignment, QC, somatic mutation calling and copy number inference
# Jingxian Clara Liu
# 2023/07/07

## clone pipeline dir
# git clone --recurse-submodules https://github.com/ding-lab/pecgs-pipeline.git
# cd pecgs-pipeline/src/compute1

## launch interactive job
# export LSF_DOCKER_VOLUMES="/storage1/fs1/dinglab/Active:/storage1/fs1/dinglab/Active /scratch1/fs1/dinglab:/scratch1/fs1/dinglab"
# export PATH="/miniconda/envs/pecgs/bin:$PATH"
# bsub -q dinglab-interactive -G compute-dinglab -Is -a 'docker(estorrs/pecgs-pipeline:0.0.1)' '/bin/bash'

RUN_LIST_ALL="/storage1/fs1/dinglab/Active/Projects/cliu/pancan_ST/PECGS/WES/WXS_request_23-05-25.tsv"
SEQ_INFO_ALL="/storage1/fs1/dinglab/Active/Projects/cliu/pancan_ST/PECGS/WES/WXS_sequencing_info_23-05-25.tsv"

## generate run commands
RUN_LIST="/storage1/fs1/dinglab/Active/Projects/cliu/pancan_ST/PECGS/WES/09.20230707_WXS4_rerun_scripts/WXS_request_23-05-25_rerun.tsv"
STORAGE_DIR="/storage1/fs1/dinglab/Active/Projects/cliu/pancan_ST/PECGS/WES/"
RUN_DIR="/scratch1/fs1/dinglab/cliu/pancan_ST/PECGS/WES/09.20230707_WXS4_rerun/"
SEQ_INFO="/storage1/fs1/dinglab/Active/Projects/cliu/pancan_ST/PECGS/WES/09.20230707_WXS4_rerun_scripts/WXS_sequencing_info_23-05-25_rerun.tsv"
SRC_DIR="/storage1/fs1/dinglab/Active/Projects/cliu/softwares/pecgs-pipeline/src/compute1/"
PIPELINE_NAME="pecgs_TN_wxs_fq"

INPUT_YAML_FILEPATH="/storage1/fs1/dinglab/Active/Projects/cliu/pancan_ST/PECGS/WES/inputs.yaml"

## [within cromwell container]
# python ${SRC_DIR}/generate_run_commands.py make-run --input-config ${INPUT_YAML_FILEPATH}  ${PIPELINE_NAME} ${RUN_LIST} ${RUN_DIR} --sequencing-info ${SEQ_INFO}
## After above step is finished:
## [outside cromwell container]
# cd ${RUN_DIR}
# bash 1.run_jobs.sh

## After above step is finished:
## [within cromwell container]
# python ${SRC_DIR}/generate_run_commands.py tidy-run ${PIPELINE_NAME} ${RUN_LIST} ${RUN_DIR} --sequencing-info ${SEQ_INFO}
## After above step is finished:
## [outside cromwell container]
# cd ${RUN_DIR}
# bash 2.tidy_run.sh

## After above step is finished:
## [within cromwell container]
# python ${SRC_DIR}/generate_run_commands.py summarize-run ${PIPELINE_NAME} ${RUN_LIST} ${RUN_DIR} --sequencing-info ${SEQ_INFO}

## After above step is finished:
## [within cromwell container]
# python ${SRC_DIR}/generate_run_commands.py move-run ${PIPELINE_NAME} ${RUN_LIST} ${RUN_DIR} --target-dir ${STORAGE_DIR}

