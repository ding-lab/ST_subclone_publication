## FUNCTION ## ---------------------------------------------------------------
getROGUE = function(obj, assay = 'Spatial', celltype_col = 'ct', sample_col = "orig.ident", span = 0.6){
    # Get expression
    expr = GetAssayData(obj[[assay]], slot = 'counts') %>% as.matrix
    #return(expr)
    # Filtering out low-abundance genes and low-quality cells
    expr <- matr.filter(expr, min.cells = 10, min.genes = 10)
    # result1: Expression entropy model
    ent.res <- SE_fun(expr)
    # result2: Get ROGUE score of the whole expression matrix
    rogue.value <- CalculateRogue(ent.res, platform = "UMI")
    # result3: Get ROGUE score of each putative cluster for each sample
    meta = obj@meta.data
    if(celltype_col %in% colnames(meta)){
        rogue.res <- rogue(expr, labels = meta[colnames(expr), celltype_col], samples = sample_col, platform = "UMI", span = span, filter= T)
    }else{
        # skip this calculation
        message('No celltype column in metadata. Skip the calculation for celltype-specific ROGUE score.')
        rogue.res = NULL
    }
    # return all 3 result values
    return(list(ent.res = ent.res, rogue.value = rogue.value, rogue.res = rogue.res))
}