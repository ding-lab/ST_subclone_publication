# Scripts related to "Focal clonal evolution in microregions" (Figure 2)

* `1_inferCNV_run.R`: Run InferCNV on Visium ST data (related to Fig. 2a-g)
    * `gencode_v32_gene_name.txt`: Generated per InferCNV instruction with reference gene annotation from `refdata-gex-GRCh38-2020-A`
    ```
    wget https://raw.githubusercontent.com/broadinstitute/infercnv/master/scripts/gtf_to_position_file.py
    python gtf_to_position_file.py --attribute_name gene_name genes.gtf gencode_v32_gene_name.txt
    ```
* `2_CalicoST_run.md`: Run CalicoST on Visiium ST data (related to Fig. 2a-g)
* `3_CNV_jaccard_similarity.R`: Calculate Jaccard similarity between two CNV profiles (related to Fig. 2a-g)
* `4_WES_CNV.sh`: Infer copy number variation from WES data (related to Fig. 2a-g)
* `5_genomic_heatmap.R`: Plot spatial CNV geonmic heatmaps (related to Fig. 2c-g)
* `6_mutation_mapping.sh`: Map WES-based somatic mutations to Visium ST (related to Fig. 2f)

