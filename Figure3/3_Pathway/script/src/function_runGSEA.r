message("conda activate clusterProfiler")
# Convert Symbol to ENTREZID 
library(clusterProfiler)
library(msigdbr)
organism = "org.Hs.eg.db"
library(organism, character.only = TRUE)
message("Using organism: ", organism)
message("Available functions: gene2entrezid, genevectorName2entrezid, genevectorName2entrezid")
gene2entrezid = function(gene, organ_database = "org.Hs.eg.db", na.rm =F){
    eg_df = bitr(gene, fromType="SYMBOL", toType="ENTREZID", OrgDb=organ_database)
    converted_id = eg_df$ENTREZID %>% setNames(eg_df$SYMBOL) 
    if(na.rm) converted_id = converted_id %>% .[!is.na(.)]
    return(converted_id)
}
genevectorName2entrezid = function(fc_vec_genelist,  organ_database = "org.Hs.eg.db", na.rm =T){
    gene_symbols = names(fc_vec_genelist)
    entrezid_vec = gene2entrezid(gene_symbols)
    gene_list_w_entrezid = fc_vec_genelist %>% setNames(entrezid_vec[gene_symbols])
    if(na.rm) gene_list_w_entrezid = gene_list_w_entrezid %>% .[!is.na(names(.))]
    return(gene_list_w_entrezid)
}

# Wrapper From DEG table to GSEA result
message("Wrapper functions: FindMarkerTable2Genelist, FindAllMarkerTable2Genelist")
# Convert from FindMarker table to gene vector.
# assume rownames already convert to column called "gene"
FindMarkerTable2Genelist = function(deg_df, gene_colname = 'gene', log2FC_colname = 'avg_log2FC'){
    # Prepare input format and remove cluster with too few DEGs
    # avg_log2FC as values, gene names as vector
    # Sort value decreasing
    deg_df %>% {setNames(.[[log2FC_colname]], .[[gene_colname]])} %>% sort(., decreasing = T)
}

# Conver from FindAllMarker table to list of gene vector 
# Assume the "cluster" column have each identity
FindAllMarkerTable2Genelist = function(deg_df, split_column = 'cluster', gene_colname = 'gene', log2FC_colname = 'avg_log2FC', n_min_deg = 5){
    # Prepare input format and remove cluster with too few DEGs
    message("Minimum number of DEG:", n_min_deg)
    deg_list= deg_df %>% 
        split(.[[split_column]]) %>% 
        # avg_log2FC as values, gene names as vector, Sort value decreasing
        map(~setNames(.x[[log2FC_colname]], .x[[gene_colname]]) %>% sort(decreasing = T)) 

    # Remove group with less than n deg
    deg_list %>% keep(~ length(.) > n_min_deg)
}   

## ------------------ GSEA ------------------
# Wrapper for Hallmark
message("Wrapper functions: FindMarkerTable2GSEAresult, FindAllMarkerTable2GSEAresult")
message("Note, default, the wrapper function uses geneset: geneset_df_HALMARK. If need replacement, please use the argument geneset_use")
message("See this page section 12.3 of how to get the dataset: http://yulab-smu.top/biomedical-knowledge-mining-book/universal-api.html")
geneset_df_HALMARK <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  #dplyr::select(gs_name, entrez_gene) %>% 
  dplyr::select(gs_name, human_gene_symbol) %>% # New, use human gene symbol
    # Remove HALLMARK
    mutate(gs_name = str_remove(gs_name, "HALLMARK_"))


ShowMsigDBHumanCategories = function(){
    message("H: hallmark gene sets
    C1: positional gene sets
    C2: curated gene sets
    C3: motif gene sets
    C4: computational gene sets
    C5: GO gene sets
    C6: oncogenic signatures
    C7: immunologic signatures")
}

LoadMSigDBHuman = function(){
    message("Loading All Human MsigDB geneset into geneset_msig_ALLHUMAN_df, geneset_msig_ALLHUMAN_use_df, geneset_msig_ALLHUMAN_list")
    geneset_msig_ALLHUMAN_df <<- msigdbr(species = "Homo sapiens") 
    geneset_msig_ALLHUMAN_list <<- geneset_msig_ALLHUMAN_df %>% split(f = .$gs_cat) %>% map(., ~dplyr::select(., gs_name, human_gene_symbol))
    ShowMsigDBHumanCategories()
    geneset_msig_ALLHUMAN_use_df <<- geneset_msig_ALLHUMAN_df %>% dplyr::select(gs_name, human_gene_symbol) 
    message("Use geneset_use = geneset_msig_ALLHUMAN_use_df in FindMarkerTable2GSEAresult")
}

FindMarkerTable2GSEAresult = function(deg_df, gene_colname = 'gene', log2FC_colname = 'avg_log2FC', geneset_use = geneset_df_HALMARK, scoreType = "std"){
    # Preflight checks ..
    gene_use = FindMarkerTable2Genelist(deg_df)
    message("Length of gene:", length(unique(deg_df[[gene_colname]])))
    #gene_vec_converted = genevectorName2entrezid(gene_use)
    # New, uses Huamn Gene symbol as reference now
    gene_vec_converted = gene_use 
    message("Length of converted:", length(gene_vec_converted))

    # first check if any gene in the geneset ID
    #geneset_overlap = intersect(geneset_use$entrez_gene, names(gene_vec_converted))
    geneset_overlap = intersect(geneset_use$human_gene_symbol, names(gene_vec_converted))
    message("Length of intersect:", length(geneset_overlap))
    if(length(geneset_overlap) == 0){message("No intersect found. Skip"); return(NULL)}
    
    # Run GSEA
    GSEA(gene_vec_converted, TERM2GENE = geneset_use, scoreType = scoreType)
}
# Run genesets from MsigDB
FindMarkerTable2GSEA_MsigDB = function(deg_df, genesets_include = c('H','C1','C2','C3','C4','C5','C6'), ...){
    ShowMsigDBHumanCategories()
    message("Using following MSigDB genesets:", toString(genesets_include))

    # Load database if not yet
    if(!exists("geneset_msig_ALLHUMAN_list" )) LoadMSigDBHuman()

    geneset_msig_list_use = geneset_msig_ALLHUMAN_list[genesets_include]
    imap(geneset_msig_list_use, function(geneset_df, cat_name){
        message("Running: ", cat_name)
        FindMarkerTable2GSEAresult(deg_df, geneset_use = geneset_df, ...)
    })
}




FindAllMarkerTable2GSEAresult = function(deg_df, split_column = 'cluster', 
    gene_colname = 'gene', log2FC_colname = 'avg_log2FC', 
    n_min_deg = 5, remove_empty_results = T,
    geneset_use = geneset_df_HALMARK,
    scoreType = "std"
    ){
    find_all_marker_gene_list = FindAllMarkerTable2Genelist(deg_df, split_column, gene_colname, log2FC_colname, n_min_deg)
    result_list = imap(find_all_marker_gene_list, function(gene_use, group_name){
        message("Running GSEA for group:", group_name)
        message("Length of gene:", length(gene_use))
        #gene_vec_converted = genevectorName2entrezid(gene_use)
        # New, uses Huamn Gene symbol as reference now
        gene_vec_converted = gene_use 

        message("Length of converted:", length(gene_vec_converted))
        # first check if any gene in the geneset IDr
        #geneset_overlap = intersect(geneset_use$entrez_gene, names(gene_vec_converted))
        geneset_overlap = intersect(geneset_use$human_gene_symbol, names(gene_vec_converted))
        message("Length of intersect:", length(geneset_overlap))
        if(length(geneset_overlap) == 0){message("No intersect found. Skip"); return(NULL)}
        print(head(gene_vec_converted))
        GSEA(gene_vec_converted, TERM2GENE = geneset_use, scoreType = scoreType)
    })

    # Clean up result
    # Could be simplify to %>% discard(.p = ~length(.)==0) if no need to print message
    if(remove_empty_results){
        message("Remove empty results")
        if(length(compact(result_list)) == 0) return(NULL)
        result_list = imap(result_list, function(gsea_result, clusterID){
            if(is.null(gsea_result)) return(NULL)
            if(nrow(gsea_result@result) == 0){
                # No GSEA result
                message("No GSEA result in ", clusterID)
                return(NULL)
            }
            message("Have GSEA result in ", clusterID)
            return(gsea_result)
        }) %>% discard(is.null)
    }  
    return(result_list)
}


# Run genesets from MsigDB
FindAllMarkerTable2GSEAresult_MsigDB = function(deg_df, genesets_include = c('H','C1','C2','C3','C4','C5','C6'), ...){
    ShowMsigDBHumanCategories()
    message("Using following MSigDB genesets:", toString(genesets_include))

    # Load database if not yet
    if(!exists("geneset_msig_ALLHUMAN_list" )) LoadMSigDBHuman()

    geneset_msig_list_use = geneset_msig_ALLHUMAN_list[genesets_include]
    imap(geneset_msig_list_use, function(geneset_df, cat_name){
        message("Running: ", cat_name)
        FindAllMarkerTable2GSEAresult(deg_df, geneset_use = geneset_df, ...)
    }) %>% discard(.p = ~length(.)==0)
}


## ---- Utility functions ----
# Make dotplot compare Tumor and TME
combineGSEAresultslist = function(gsea_result_list){
    gsea_result_list %>% imap(~as.data.frame(.x) %>% rownames_to_column('geneset') %>% mutate(group = .y) 
    )  %>% bind_rows() 
}

ReorderByHCluster = function(data_use, ident_column, groupby_column, values_column){
    # 1. convert to wide
    value_mtx = data_use %>% 
        dplyr::select(all_of(c(ident_column,groupby_column, values_column))) %>%
        pivot_wider(id_cols = all_of(ident_column), names_from = all_of(groupby_column), values_from = all_of(values_column), values_fill=0) %>% 
        column_to_rownames(ident_column) %>% 
        as.matrix() %>% 
        t()

    # safety, check if more than 2 rows. if not return as is
    if(nrow(data_use) < 2){
        message("Less than 2 row. Cannot cluster. return data as is")
        return(data_use)
    }
    # 2. cluster using hclust
    hcluster = hclust(dist(value_mtx))
    
    # 3. get the order
    label_ordered = hcluster$labels[hcluster$order]
    
    # 4. reorder the data
    data_use = data_use %>% mutate( {{ groupby_column}} := factor(.data[[groupby_column]], levels = label_ordered))
    message("Reordered groupby_column: ", groupby_column, " by hclust")
    # 5. return the data
    return(data_use)
}

# Clean geneset name
CleanGenesetName = function(geneset_names, terms_to_remove = c("HALLMARK"), convert_to_space = T, convert_to_title = T){
    search_pattern = paste(str_c(terms_to_remove,"_"), collapse = "|")
    cleaned_string = str_remove_all(geneset_names, search_pattern) 
    if(convert_to_space) cleaned_string = str_replace_all(cleaned_string, "_", " ")
    if(convert_to_title) cleaned_string = str_to_title(cleaned_string)
    return(cleaned_string)
}

# For A single GSEA result --------------------------------
GetGSEAgenes = function(gsea_result){
    gsea_result %>% as.data.frame %>% 
        dplyr::select(ID, core_enrichment) %>% 
        separate_rows(core_enrichment, sep = '/') 
}

GetTopGSEAgeneList = function(gsea_result, top_n = NULL){
    gsea_df = GetGSEAgenes(gsea_result)
    gsea_list = split(x = gsea_df$core_enrichment, f = gsea_df$ID) 
    if(is.null(top_n)){
        message("No top_n specified, return all")
    }else{
        message("Return top_n:", top_n)
        gsea_list %>% map(.f = ~head(.x, top_n))
    }
    return(gsea_list)
}

# For A LIST OF GSEA result --------------------------------
# Get summary of values based on the values_columns (default: NES, qvalue, rank)
SummarizeGSEAlist = function(gsea_use, value_columns = c('NES','qvalue', 'rank'), exclude_idents=NULL){
    summary_df = imap(gsea_use, function(gsea_result, clone_id){
        gsea_result %>% 
            as.data.frame %>% 
            dplyr::select(ID, all_of(value_columns)) %>% 
            mutate(clone_id = clone_id)
    }) %>% bind_rows() %>% remove_rownames()

    # Filter out and summarize NES
    summary_df %>% filter(!clone_id %in% exclude_idents) %>% 
        group_by(ID) %>%
        summarize(across(value_columns, mean)) %>% 
        arrange(desc(across(value_columns)))%>%  # https://rlang.r-lib.org/reference/topic-data-mask-programming.html
        as.data.frame
}


## ---- Plot functions ----
MakeGeneSetDotplot = function(gsearesult_list, plt_title = "GSEA result", subtitle = "Geneset", reorder_geneset = T, reorder_group = T, clean_id = T, ...){
    combined_df = combineGSEAresultslist(gsearesult_list)
    # Clean geneset name
    if(clean_id) combined_df = combined_df %>% mutate(geneset = CleanGenesetName(geneset))
    # Make dotplot x = group, y = geneset, color = enrichmentScore, size = -log10(qvalue)
    # Reorder geneset/group by hclust
    if(reorder_geneset) combined_df = ReorderByHCluster(data_use = combined_df, 'group', 'geneset', 'NES') 
    if(reorder_group) combined_df = ReorderByHCluster(data_use = combined_df, "geneset", 'group', 'NES')
        
    ggplot(combined_df, aes(x = group, y = geneset, color = NES, size = -log10(qvalue))) + geom_point() + 
        scale_color_gradient2(low = "#2222EE", mid = 'gray80', high = "#FF4444") + 
        theme_bw() + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
        #facet_wrap(~group, scales = 'free_y') + 
        labs(title = plt_title, subtitle = subtitle, color = "Normalized\nEnrichment\nScore", size = '-log10(qvalue)') + 
        coord_fixed()
        #theme(legend.position = 'none')
}
