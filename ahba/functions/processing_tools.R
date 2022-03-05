#' Aggregate AHBA samples using atlas mappings
#'
#' @param donorInfo (character) Vector with elements 1. the donor name, 2. the path to the preprocessed data, 3. the path to the atlas mappings
#' @param atlas (character) String indicating which atlas mapping is used
#'
#' @return
aggregateSamples <- function(donorInfo, atlas = "AHBA"){
  
  donorID <- donorInfo[[1]]
  dataPath <- donorInfo[[2]]
  atlasPath <- donorInfo[[3]]
  
  #Read in preprocessed expression data  
  dfMicroarrayExpression <- suppressMessages(read_csv(dataPath))
  
  #Read in atlas mappings
  dfAtlasMappings <- suppressMessages(read_csv(atlasPath))
  
  if(atlas == "AALsubcortical"){
    
    dfAtlasMappings <- dfAtlasMappings %>% 
      select(SampleID, 
             AtlasStructure,
             AtlasLabel) %>% 
      filter(AtlasLabel != 0)
    
  } else {
    
    dfAtlasMappings <- dfAtlasMappings %>% 
      select(SampleID,
             AtlasStructure = structure_name, 
             AtlasLabel = structure_id)
    
  }
  
  #Aggregate samples per structure
  dfMicroarrayExpression_wLabels <- dfMicroarrayExpression %>% 
    gather(key = SampleID, value = Expression, -1) %>% 
    left_join(dfAtlasMappings, by = "SampleID") %>% 
    mutate(Expression = 2^Expression) %>% 
    group_by(Gene, AtlasStructure) %>% 
    summarise(StructExpression = mean(Expression, na.rm = T)) %>% 
    ungroup() %>% 
    mutate(StructExpression = log2(StructExpression)) %>% 
    spread(key = AtlasStructure, value = StructExpression, -1)
  
  
  #Order the columns according to the atlas labels
  structOrder <- dfAtlasMappings %>% 
    select(AtlasStructure, AtlasLabel) %>% 
    unique() %>% 
    arrange(AtlasLabel) %>% 
    pull(AtlasStructure)
  
  structOrder <- structOrder[structOrder %in% colnames(dfMicroarrayExpression_wLabels[,-1])]
  indOrder <- c(1, match(structOrder, colnames(dfMicroarrayExpression_wLabels)))
  
  dfMicroarrayExpression_wLabels <- dfMicroarrayExpression_wLabels[, indOrder]
  
  return(dfMicroarrayExpression_wLabels)
  
}


#' Combine donor gene data into one
#'
#' @param data (list)
#' @param verbose (logical)
#'
#' @return
combineDonorGenes <- function(data, verbose = TRUE){
  
  if(verbose) {message("Combining donor genes...")}
  
  genesShared <- map(data, function(x){x[["GeneExpression"]][["Gene"]]}) %>% 
    reduce(intersect) %>% 
    sort()
  
  dfGeneExpression <- map_dfc(data, 
                              function(x){
                                x[["GeneExpression"]] %>% 
                                  filter(Gene %in% genesShared) %>% 
                                  arrange(Gene) %>% 
                                  select(-Gene)
                              }) %>% 
    mutate(Gene = genesShared)
  
  #Combine sample information into one data frame
  dfSampleInfo <- map_dfr(data, 
                          function(x){x[["SampleInfo"]]})  
  
  #Store donor names in a vector
  donors <- map_chr(data,
                    function(x){x[["Donor"]]})
  names(donors) <- NULL
  
  return(list(Donor = donors,
              GeneExpression = dfGeneExpression,
              SampleInfo = dfSampleInfo))
}


#' Combine donor probe data into one
#'
#' @param data (list)
#' @param verbose (logical)
#'
#' @return
combineDonorProbes <- function(data, verbose = TRUE){
  
  if(verbose) {message("Combining donor probes...")}
  
  #Extract probe information (same for all donors)
  dfProbeInfo <- data[[1]][["ProbeInfo"]]
  
  #Combine probe expression values into one data frame
  dfProbeExpression <- map_dfc(data,
                               function(x){select(x[["ProbeExpression"]], -ProbeID)}) %>% 
    mutate(ProbeID = dfProbeInfo$ProbeID)
  
  dfIntensityFilter <- map_dfc(data,
                               function(x){select(x[["IntensityFilter"]], -ProbeID)}) %>% 
    mutate(ProbeID = dfProbeInfo$ProbeID)
  
  #Combine sample information into one data frame
  dfSampleInfo <- map_dfr(data, 
                          function(x){x[["SampleInfo"]]})
  
  #Store donor paths in a vector
  paths <- map_chr(data, 
                   function(x){x[["DonorPath"]]})
  
  #Store donor names in a vector
  donors <- map_chr(data,
                    function(x){x[["Donor"]]})
  names(donors) <- NULL
  
  return(list(Donor = donors,
              DonorPath = paths,
              ProbeInfo = dfProbeInfo,
              ProbeExpression = dfProbeExpression,
              SampleInfo = dfSampleInfo,
              IntensityFilter = dfIntensityFilter))
}


#' Extract samples from the left hemisphere
#'
#' @param data (list)
#'
#' @return
extractLeftHemisphere <- function(data){
  
  data[["SampleInfo"]] <- data[["SampleInfo"]] %>% 
    filter(mni_x < 0)
  
  indLeft <- colnames(data[["ProbeExpression"]]) %in% c("ProbeID", data[["SampleInfo"]]$SampleID)
  data[["ProbeExpression"]] <- data[["ProbeExpression"]][,indLeft]
  
  return(data)
  
}


#' Filter the AHBA microarray probes
#'
#' @param data (list)
#' @param entrezFiltering (logical) Flag for entrez ID filtering
#' @param intensityFiltering (logical) Flag for intensity filtering
#' @param intensityThreshold (numeric) Threshold to use for intensity filtering
#' @param verbose (logical)
#'
#' @return
filterDonorData <- function(data, entrezFilter = TRUE, intensityFilter = TRUE, intensityThreshold = 0.5, verbose = TRUE){
  
  if(verbose) {message("Filtering microarray probes...")}
  
  if(entrezFilter){
    #Define and apply entrezID filter
    if(verbose) {message("Applying entrez ID filter...")}
    entrezFilterPass <- !is.na(data[["ProbeInfo"]]$EntrezID)
    data[["ProbeInfo"]] <- data[["ProbeInfo"]][entrezFilterPass,]
    data[["ProbeExpression"]] <- data[["ProbeExpression"]][entrezFilterPass,]
    data[["IntensityFilter"]] <- data[["IntensityFilter"]][entrezFilterPass,]
  }
  
  if(intensityFilter){
    
    if(verbose) {message("Applying intensity filter...")}
    
    #Convert intensity filter to matrix
    matIntensityFilter <- data[["IntensityFilter"]] %>% 
      select(-ProbeID) %>% 
      as.matrix()
    
    #Index for probes that pass the threshold
    intensityFilterPass <- rowSums(matIntensityFilter)/ncol(matIntensityFilter) > intensityThreshold
    
    #Apply the filter
    data[["ProbeInfo"]] <- data[["ProbeInfo"]][intensityFilterPass,]
    data[["ProbeExpression"]] <- data[["ProbeExpression"]][intensityFilterPass,]
    data[["IntensityFilter"]] <- data[["IntensityFilter"]][intensityFilterPass,]
  }
  
  return(data)
  
}


#' Import the AHBA microarray expression data for a given donor
#'
#' @param path (character) Path to donor data directory
#' @param donor (character) Name of donor
#' @param verbose (logical)
#'
#' @return 
importDonorData <- function(path, donor, verbose = TRUE){
  
  if(verbose) {message(str_c("Importing data for donor ", donor, "..."))}
  
  #Microarray probe (rows) expression across all samples (cols)
  dfProbeExpression <- suppressMessages(read_csv(str_c(path, "MicroarrayExpression.csv"), col_names = F)) %>% 
    rename(ProbeID = X1)
  
  #Microarray probe information
  dfProbeInfo <- suppressMessages(read_csv(str_c(path, "Probes.csv"))) %>% 
    select(ProbeID = probe_id,
           Gene = gene_symbol,
           EntrezID = entrez_id)
  
  #Sample information
  dfSampleInfo <- suppressMessages(read_csv(str_c(path, "SampleAnnot.csv"))) %>% 
    mutate(SampleID = str_c(structure_id, slab_num, well_id, sep = "-"),
           Donor = donor)
  
  dfIntensityFilter <- suppressMessages(read_csv(str_c(path, "PACall.csv"), col_names = F)) %>% 
    rename(ProbeID = X1)
  
  #Label expression data columns with sample IDs
  colnames(dfProbeExpression)[2:ncol(dfProbeExpression)] <- dfSampleInfo$SampleID
  colnames(dfIntensityFilter)[2:ncol(dfIntensityFilter)] <- dfSampleInfo$SampleID
  
  return(list(Donor = donor,
              DonorPath = path,
              ProbeExpression = dfProbeExpression,
              ProbeInfo = dfProbeInfo,
              SampleInfo = dfSampleInfo,
              IntensityFilter = dfIntensityFilter))
}


#' Process the AHBA data
#'
#' @param data 
#' @param version 
#' @param verbose 
#'
#' @return
processingPipeline <- function(data, version = 1, verbose = TRUE){
  if (version == 1){
    message("Running pipeline version 1")
    data %>% 
      map(filterDonorData,
          entrezFilter = TRUE,
          intensityFilter = TRUE,
          intensityThreshold = 0.5,
          verbose = verbose) %>% 
      map(selectMicroarrayProbes,
          method = "Beauchamp",
          verbose = verbose) %>% 
      combineDonorGenes(verbose = verbose)
  } else if (version == 2){
    message("Running pipeline version 2")
    data %>% 
      combineDonorProbes(verbose = verbose) %>% 
      filterDonorData(entrezFilter = TRUE,
                      intensityFilter = TRUE,
                      intensityThreshold = 0.5,
                      verbose = verbose) %>% 
      selectMicroarrayProbes(method = "Beauchamp",
                             verbose = verbose)
  } else {
    stop("Version not valid")
  }
}


#' Consolidate multiple microarray probes for each gene
#'
#' @param data (list)
#' @param method (character) One of "Beauchamp" or "Myers", indicating how to consolidate probes per gene
#' @param verbose (logical)
#'
#' @return
selectMicroarrayProbes <- function(data, method = "Beauchamp", verbose = TRUE){
  
  if(verbose){message("Consolidating multiple probes per gene...")}
  
  if (method == "Beauchamp") {
    
    #Average expression of multiple probes per gene (on linear scale)
    expr <- data[["ProbeExpression"]] %>% 
      inner_join(data[["ProbeInfo"]],
                 by = "ProbeID") %>% 
      select(-EntrezID, -ProbeID) %>%
      mutate_at(.vars = vars(-Gene), .funs = function(x){2^x}) %>% 
      group_by(Gene) %>% 
      summarise_all(.funs = mean) %>% 
      ungroup() %>%
      mutate_at(.vars = vars(-Gene), .funs = function(x){log2(x)}) %>% 
      arrange(Gene)
    
  } else if (method == "Myers"){
    
    #Number of probes per gene
    probeCount <- data[["ProbeInfo"]] %>% 
      group_by(Gene) %>% 
      summarise(nProbes = n())
    
    #Genes with 1 probe
    genes_1 <- probeCount %>% 
      filter(nProbes == 1) %>% 
      pull(Gene)
    
    #Select unique probe
    probes_1 <- data[["ProbeInfo"]] %>% 
      filter(Gene %in% genes_1) %>% 
      pull(ProbeID)
    
    #Genes with 2 probes
    genes_2 <- probeCount %>% 
      filter(nProbes == 2) %>% 
      pull(Gene)
    
    df_2 <- data[["ProbeInfo"]] %>% 
      filter(Gene %in% genes_2) %>% 
      select(Gene, ProbeID)
    
    #Select probe per gene based on maximum mean intensity across all samples
    probes_2 <- data[["ProbeExpression"]] %>% 
      filter(ProbeID %in% df_2$ProbeID) %>% 
      column_to_rownames("ProbeID") %>% 
      as.matrix() %>% 
      rowMeans() %>% 
      enframe(name = "ProbeID",
              value = "MeanIntensity") %>% 
      mutate(ProbeID = as.numeric(ProbeID)) %>% 
      inner_join(df_2, by = "ProbeID") %>% 
      group_by(Gene) %>% 
      filter(MeanIntensity == max(MeanIntensity)) %>% 
      ungroup() %>% 
      pull(ProbeID)
    
    #Genes with > 2 probes
    genes_many <- probeCount %>% 
      filter(nProbes > 2) %>% 
      pull(Gene) 
    
    df_many <- data[["ProbeInfo"]] %>% 
      filter(Gene %in% genes_many) %>% 
      select(Gene, ProbeID)
    
    #Select probe per gene based on maximal average correlation to all other probes for that gene
    probes_many <- data[["ProbeExpression"]] %>% 
      inner_join(df_many, by = "ProbeID") %>% 
      group_by(Gene) %>% 
      nest() %>% 
      mutate(MeanCorrelation = map(data, 
                                   function(x){
                                     x %>% 
                                       column_to_rownames("ProbeID") %>% 
                                       as.matrix() %>% 
                                       t() %>% 
                                       cor() %>% 
                                       rowMeans()
                                   }),
             ProbeID = map(MeanCorrelation, ~ names(.x))) %>% 
      select(Gene, MeanCorrelation, ProbeID) %>% 
      unnest(cols = c("MeanCorrelation", "ProbeID")) %>% 
      group_by(Gene) %>% 
      filter(MeanCorrelation == max(MeanCorrelation)) %>% 
      ungroup() %>% 
      pull(ProbeID)
    
    #Combine probe selections 
    probes <- c(probes_1, probes_2, probes_many)
    
    #Extract selected probes
    indProbes <- data[["ProbeInfo"]]$ProbeID %in% probes
    expr <- inner_join(data[["ProbeInfo"]][indProbes,],
                       data[["ProbeExpression"]][indProbes,],
                       by = "ProbeID") %>% 
      select(-ProbeID, -EntrezID) %>% 
      arrange(Gene)
    
  } else {
    stop()
  }
  
  return(list(Donor = data[["Donor"]],
              GeneExpression = expr,
              SampleInfo = data[["SampleInfo"]]))
  
}