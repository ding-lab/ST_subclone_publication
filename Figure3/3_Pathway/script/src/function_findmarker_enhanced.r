library(tidyverse)
library(Seurat)

message("Functions includes: FindMarkersFormated, FindMarkersEachVsRef, FindMarkersEachVsRefComplete")
message("FindMarkers with better format: FindMarkersFormated")
FindMarkersFormated = function(obj, group.by = 'seurat_clusters', ident.1=NULL, ident.2=NULL, ...){
    message("group.by use: ", group.by)
    Idents(obj) = group.by
    message("ident.1 names as 'cluster' to be same as FindAllMarkers")
    FindMarkers(obj,ident.1 = ident.1, ident.2=ident.2,  ...) %>% 
        rownames_to_column('gene') %>%
        mutate(cluster = as.character(ident.1),
            ident2 = ident.2,
            group_by = group.by
        ) 
}

# Compare every other ident in the group.by to ref
message("Compare every other ident in the group.by to ref: FindMarkersEachVsRef")
FindMarkersEachVsRef = function(obj, group.by='seurat_clusters', ident_ref, ...){
    message('Ref. ident is: ', ident_ref)
    remain_idents = setdiff(obj@meta.data[[group.by]], ident_ref)
    message('Remain idents are: ', toString(remain_idents))
    map(remain_idents, function(ident_1_use){
        message('Ident1 use: ', ident_1_use)
        FindMarkersFormated(obj, group.by = group.by, ident.1 = ident_1_use, ident.2 = ident_ref, ...)
    }) %>% bind_rows()
}

# Compare every other ident in the group.by to ref
# Also include Ref vs NotRef
message("Compare every other ident in the group.by to ref, also include Ref vs NotRef: FindMarkersEachVsRefComplete")
FindMarkersEachVsRefComplete = function(obj, group.by='seurat_clusters', ident_ref, ...){
    message('Ref. ident is: ', ident_ref)
    # 1. Ref vs Not Ref
    deg_ref_notRef_df = FindMarkersFormated(ST_use, group.by = group.by, ident.1 = ident_ref) %>% 
        mutate(ident2 = str_glue("Not_{ident_ref}"))
    
    # 2. Each other ident vs Ref
    remain_idents = setdiff(obj@meta.data[[group.by]], ident_ref)
    message('Remain idents are: ', toString(remain_idents))
    deg_each_vs_ref_df = map(remain_idents, function(ident_1_use){
        message('Ident1 use: ', ident_1_use)
        FindMarkersFormated(obj, group.by = group.by, ident.1 = ident_1_use, ident.2 = ident_ref, ...)
    }) %>% bind_rows()

    # 3. Combine
    return(bind_rows(deg_ref_notRef_df, deg_each_vs_ref_df))
}

# Convert existing find marker result to be same as FindMarkersFormated
# EXTERNAL
# Change ident1 to cluster, and Gene/Features to gene
Convert2FindMarkersFormated = function(deg_df, ident1_colname = 'ident.1', gene_colname = 'Gene', ...){
    ## ---------- IDENT1 ------------- ##
    message("This function convert ident.1 names as 'cluster' to be same as FindAllMarkers")
    message("This function convert Gene/Features names as 'gene' to be same as FindAllMarkers")
    deg_df %>% Convert2FindMarkers_(
        from_colname = ident1_colname, 
        to_colname = 'cluster', 
        guess_pattern = '^[I,i]dent.?1$', ...) %>% 
            Convert2FindMarkers_(
                from_colname = gene_colname, 
                to_colname = 'gene', 
                guess_pattern = '^[G,g]ene[s]?$|^[F,f]eature[s]?$|^[S,s]ymbol[s]?',
                ...)
}

# Internal function
Convert2FindMarkers_ = function(deg_df, from_colname=NULL, to_colname, guess_pattern = '^[I,i]dent.?1$'){
    ## ---------- IDENT1 ------------- ##
    message("This function convert ", from_colname, " names to ", to_colname,  " to be same as FindAllMarkers")
    if(to_colname %in% names(deg_df)){
        message('Already have target column nothing to convert')
        return(deg_df)
    }
    if(from_colname %in% names(deg_df)){ # name provided. Directly
        message("Convert ", from_colname, " to ", to_colname)
        deg_df = deg_df %>% 
            dplyr::rename({{to_colname}} := all_of(from_colname))
        return(deg_df)
    }
    # Name not provided or not in the data frame ; guess
    # Guess the ident1 name variant: possibles: ident.1, ident1, Ident1, Ident.1, IDENT.1 ...
    df_column_names = names(deg_df)
    from_colname = df_column_names[df_column_names %>% tolower() %>% str_detect(guess_pattern)]
    if(length(from_colname) != 1){
    message("column name guessing failed.. please provide the column name")
    return(deg_df)
    }
    # rename
    deg_df = deg_df %>% dplyr::rename(
        {{to_colname}} := !!sym(from_colname)
    )
    message("Renamed ", from_colname, " to ", to_colname)
    return(deg_df)
}
