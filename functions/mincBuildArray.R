#' Build a MINC array from structure-wise data
#' 
#' @description 
#' Map a vector of structure-wise values to voxels in an image.
#' 
#' @details 
#' This function matches the values defined for every region of 
#' interest (ROI) in an atlas to the voxels that lie within those 
#' ROIs in an image. This is done either by matching ROI labels 
#' directly to an atlas, or by matching ROI definitions to atlas labels 
#' using an intermediate dictionary.
#' 
#' @param values (numeric vector) A named vector of some quantity for 
#' each region in an atlas. Names can be atlas labels or atlas label 
#' definitions.
#' @param labels (mincSingleDim) An array containing atlas labels at 
#' every voxel.
#' @param defs (data.frame) A data frame containing columns "Structure" 
#' and "Label" indicating the mapping between atlas labels and 
#' definitions. Can be left as `NULL` if `names(values)` contains atlas 
#' labels, but must be defined if `names(values)` contains atlas 
#' definitions.
#' @param values.names (character vector) Optional vector of names to 
#' use for `names(values)`. If not `NULL`, will override whatever names 
#' values already has.
#' 
#' @return (mincSingleDim) An array with structure-wise values mapped 
#' to voxels.
mincBuildArray <- function(values, labels, defs = NULL, values.names = NULL){
  
  #If values.names is passed, assign those names to values
  if(!is.null(values.names)){
    names(values) <- values.names
  } 
  
  #Check: values must be numeric
  if(!is.numeric(values)){
    stop(paste("Argument `values` must be a numeric vector. Got class:", 
               class(values)))
  }
  
  #Check: values must be a named vector
  if(is.null(names(values))){
    stop(paste("Argument `values` must be a named vector with names",
               "corresponding to atlas labels or definitions. Names",
               "can be assigned prior to function call or passed using",
               "argument `values.names`"))
  }
  
  #Control flow for names being labels or label definitions
  numericNames <- suppressWarnings(as.numeric(names(values)))
  if (sum(is.na(numericNames)) == 0) {
    
    #If some values names are not in the atlas, they will not be mapped
    #to any voxels
    notInAtlas <- names(values)[!(as.numeric(names(values)) %in% labels)]
    if(length(notInAtlas) != 0){
      warning(paste("The entries in `values` with the following names have no",
                    "corresponding atlas labels in `labels`: ",
                    paste(notInAtlas, collapse = ", "),
                    ". They will not be mapped to any voxels in the atlas."))
    }
    
    #Identify atlas labels that aren't in values. 
    #These will have to be set to 0.
    notInValues <- unique(labels)[!(unique(labels) %in% as.numeric(names(values)))]
    
    #Remove empty voxels, which aren't in the atlas anyway
    notInValues <- notInValues[notInValues != 0] 
    
    if(length(notInValues) != 0){
      warning(paste("The following atlas labels were not associated with any",
                    "entries in `values`: ", 
                    paste(notInValues, collapse = ", "),
                    ". They will be mapped to zero."))
    }
    
    #Assign a value of 0 to those atlas labels not in values
    valuesMissing <- rep(0, length(notInValues))
    names(valuesMissing) <- notInValues
    
    #Add the 0 values
    values <- c(values, valuesMissing)
    
    #Map the values onto the voxels using atlas labels
    out <- values[match(labels, as.numeric(names(values)))]
    
  } else {
    
    #When the names of values are atlas definitions, 
    #we need a mapping between definitions and labels
    if(is.null(defs)){
      stop(paste("Argument `defs` must be supplied when argument `values` has",
                 "names corresponding to atlas definitions instead of labels",
                 "(checked by coercing names to numeric)"))
    }
    
    #Check: defs must be a data frame
    if(!is.data.frame(defs)){
      stop(paste("Argument `defs` must be a data frame. Got class:", 
                 class(defs)))
    }
    
    #Check: defs must have columns Structure and Label
    if(any(!(c("Structure", "Label") %in% colnames(defs)))){
      stop(paste("Argument `defs` should have columns `Structure` and `Label`",
                 "containing atlas definitions and labels, respectively"))
    }
    
    #If some values names are not in the atlas, they will not be mapped
    #to any voxels
    notInAtlas <- names(values)[!(names(values) %in% defs$Structure)]
    if(length(notInAtlas) != 0){
      warning(paste("The entries in `values` with the following names have no",
                    "corresponding entries in `defs`: ",
                    paste(notInAtlas, collapse = ", "),
                    ". They will not be mapped to any voxels in the atlas."))
    }
    
    #Identify atlas regions not found in values. These must be mapped to 0.
    notInValues <- defs$Structure[!(defs$Structure %in% names(values))]
    
    if(length(notInValues) != 0){
      warning(paste("The following atlas regions were not associated with any",
                    "entries in `values`: ", 
                    paste(notInValues, collapse = ", "),
                    ". They will be mapped to zero."))
    }
    
    #Assign a value of 0 to those missing structures
    valuesMissing <- rep(0, length(notInValues))
    names(valuesMissing) <- notInValues
    
    #Add the 0 value to the vector of interest
    values <- c(values, valuesMissing)
    
    #Match the vector of values to the atlas labels
    defs$Values <- values[defs$Structure] 
    
    #Map the values onto the voxels using atlas labels
    out <- defs$Values[match(labels, defs$Label)]
    
  }
  
  #Set empty voxels to 0
  out[is.na(out)] <- 0
  
  #Copy attributes of labels array
  attributes(out) <- attributes(labels)
  
  return(out)
}