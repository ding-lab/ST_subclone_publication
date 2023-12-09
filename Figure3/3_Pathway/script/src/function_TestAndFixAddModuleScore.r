### Fix Seurat object for AddModuleScore
# Detail: Fix the Error in `cut_number()`:  ! Insufficient data values to produce 24 bins.
# issue by removing features with less than 2 cells 
TestAndFixSeuratForAddModuleScore = function(obj){
    # 1. Create a modified AddModuleScore function that returns NULL if there is an error
    AddModuleScorePossibly = possibly(
        .f = AddModuleScore,
        otherwise = NULL)

    # 2. Run test the function 
    test_result = obj %>% 
        AddModuleScorePossibly(., features = rownames(.)[1:3])
    
    # 3. Check if need to fix the object
    if(is.null(test_result)){
        message("Got error running AddModuleScore! might have feature with <= 1 cells. Attempt to use RemoveZeroCellFeatures_ to fix the object")
        obj_fixed = RemoveZeroCellFeatures_(obj)
        test_result_final = obj_fixed %>% 
            AddModuleScorePossibly(., features = rownames(.)[1:3])
        if(is.null(test_result_final)){
            message("Unable to fix the object, please check the seurat object")
            return(obj)
        }
        message("Successfully fixed the object")
        return(obj_fixed)
    }
    message("Nothing to fix!")
    return(obj)
        
}

# Try subset and do again V2 works!
RemoveZeroCellFeatures_ = function(obj, assays){
    if(missing(assays)) assays = Assays(obj)
    message("Removing zero cell features in assays:", toString(assays))
    obj_new = obj
    for(assay_use in assays){
        message("Processing assay:", assay_use)
        message("N features before filtering:", nrow(obj_new@assays[[assay_use]]))
        # Find features to keep
        features_more_cell = st_obj_use@assays[[assay_use]] %>% GetAssayData %>% rowSums() %>% .[. >1 ] %>% names
        obj_new[[assay_use]] = 
            subset(obj_new[[assay_use]], features = features_more_cell)
                
        #st_obj_use@assays[[assay_use]][features_more_cell,]
        message("N features after filtering:", nrow(obj_new@assays[[assay_use]]))
    }
    return(obj_new)
}
