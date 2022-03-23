library(tidyverse)

#' Combine donor gene data into one
#'
#' @param data (list) A list containing elements:
#'              Donor (character) 
#'                    Donor name
#'                GeneExpression (tibble)
#'                    Tibble containing gene-by-sample expression values
#'                SampleInfo (tibble) 
#'                    Tibble containing metadata for each sample
#' @param verbose (logical) Verbosity option
#'
#' @return (list) A list containing three elements:
#'                Donor (character) 
#'                    The names of all donors
#'                GeneExpression (tibble)
#'                    Tibble containing gene-by-sample expression values
#'                    for all donors
#'                SampleInfo (tibble) 
#'                    Tibble containing metadata for each sample for 
#'                    all donors
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
#' @param data (list) A list with elements:
#'                Donor (character) 
#'                    Donor name
#'                DonorPath (character)
#'                    Path to donor data dir 
#'                ProbeExpression (tibble) 
#'                    Tibble containing probe-by-sample expression
#'                ProbeInfo (tibble) 
#'                    Tibble containing metadata for each probe
#'                SampleInfo (tibble) 
#'                    Tibble containing metadata for each sample
#'                IntensityFilter (tibble)
#'                    Tibble containing probe-by-sample binary filter/mask
#' @param verbose (logical)
#'
#' @return (list) A list containing six elements:
#'                Donor (character) 
#'                    The names of all donors
#'                DonorPath (character)
#'                    Paths to each donor data dir 
#'                ProbeExpression (tibble) 
#'                    Tibble containing probe-by-sample expression
#'                ProbeInfo (tibble) 
#'                    Tibble containing metadata for each probe for all donors
#'                SampleInfo (tibble) 
#'                    Tibble containing metadata for each sample for all donors
#'                IntensityFilter (tibble)
#'                    Tibble containing probe-by-sample binary filter/mask for 
#'                    all donors
combineDonorProbes <- function(data, verbose = TRUE){
  
  if(verbose) {message("Combining donor probes...")}
  
  #Extract probe information (same for all donors)
  dfProbeInfo <- data[[1]][["ProbeInfo"]]
  
  #Combine probe expression values into one data frame
  dfProbeExpression <- map_dfc(data,
                               function(x){
                                 select(x[["ProbeExpression"]], -ProbeID)
                               }) %>% 
    mutate(ProbeID = dfProbeInfo$ProbeID)
  
  dfIntensityFilter <- map_dfc(data,
                               function(x){
                                 select(x[["IntensityFilter"]], -ProbeID)
                               }) %>% 
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


#' Filter the AHBA microarray probes
#'
#' @param data (list) A list with elements:
#'                Donor (character) 
#'                    Donor name
#'                DonorPath (character)
#'                    Path to donor data dir 
#'                ProbeExpression (tibble) 
#'                    Tibble containing probe-by-sample expression
#'                ProbeInfo (tibble) 
#'                    Tibble containing metadata for each probe
#'                SampleInfo (tibble) 
#'                    Tibble containing metadata for each sample
#'                IntensityFilter (tibble)
#'                    Tibble containing probe-by-sample binary filter/mask
#' @param entrezFiltering (logical) Entrez ID filtering
#' @param intensityFiltering (logical) Intensity filtering
#' @param intensityThreshold (numeric) Threshold to use for intensity filtering
#' @param verbose (logical) Verbosity option
#'
#' @return A list with the same components as `data`, but with filters applied.
filterDonorData <- function(data, entrezFilter = TRUE, intensityFilter = TRUE, 
                            intensityThreshold = 0.5, verbose = TRUE){
  
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
#' @param verbose (logical) Verbosity option
#'
#' @return (list) A list with elements:
#'                Donor (character) 
#'                    Donor name
#'                DonorPath (character)
#'                    Path to donor data dir 
#'                ProbeExpression (tibble) 
#'                    Tibble containing probe-by-sample expression
#'                ProbeInfo (tibble) 
#'                    Tibble containing metadata for each probe
#'                SampleInfo (tibble) 
#'                    Tibble containing metadata for each sample
#'                IntensityFilter (tibble)
#'                    Tibble containing probe-by-sample binary filter/mask
importDonorData <- function(path, donor, verbose = TRUE){
  
  if(verbose) {message(str_c("Importing data for donor ", donor, "..."))}
  
  #Microarray probe (rows) expression across all samples (cols)
  pathProbeExpr <- str_c(path, "MicroarrayExpression.csv")
  dfProbeExpression <- suppressMessages(read_csv(pathProbeExpr, 
                                                 col_names = F)) %>% 
    rename(ProbeID = X1)
  
  #Microarray probe information
  pathProbeInfo <- str_c(path, "Probes.csv")
  dfProbeInfo <- suppressMessages(read_csv(pathProbeInfo)) %>% 
    select(ProbeID = probe_id,
           Gene = gene_symbol,
           EntrezID = entrez_id)
  
  #Sample information
  pathSampleInfo <- str_c(path, "SampleAnnot.csv")
  dfSampleInfo <- suppressMessages(read_csv(pathSampleInfo)) %>% 
    mutate(SampleID = str_c(structure_id, slab_num, well_id, sep = "-"),
           Donor = donor)
  
  pathIntensityFilter <- str_c(path, "PACall.csv")
  dfIntensityFilter <- suppressMessages(read_csv(pathIntensityFilter,
                                                 col_names = F)) %>% 
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
#' @param data (list) A list with elements:
#'        Donor (character) 
#'            Donor name
#'        DonorPath (character)
#'            Path to donor data dir 
#'        ProbeExpression (tibble) 
#'            Tibble containing probe-by-sample expression
#'        ProbeInfo (tibble) 
#'            Tibble containing metadata for each probe
#'        SampleInfo (tibble) 
#'            Tibble containing metadata for each sample
#'        IntensityFilter (tibble)
#'            Tibble containing probe-by-sample binary filter/mask
#' @param version (numeric scalar) Pipeline version (1, 2)
#'        Version 1
#'            1. Filter data for each donor
#'            2. Average multiple probes per gene for each donor
#'            3. Combine donor samples into one matrix
#'        Version 2
#'        1. Combine probes from all donors into one matrix
#'        2. Filter data for all donors together
#'        3. Select one probe per gene using Myers' 3-step method         
#' @param verbose (logical scalar)
#'
#' @return (list) A list containing three elements:
#'         Donor (character) 
#'             The names of all donors
#'         GeneExpression (tibble)
#'             Tibble containing gene-by-sample expression values
#'             for all donors
#'         SampleInfo (tibble) 
#'             Tibble containing metadata for each sample for 
#'             all donors
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
#' @param data (list) A list containing elements:
#'                Donor (character) 
#'                    Donor name
#'                ProbeExpression (tibble) 
#'                    Tibble containing probe-by-sample expression
#'                ProbeInfo (tibble) 
#'                    Tibble containing metadata for each probe
#'                SampleInfo (tibble) 
#'                    Tibble containing metadata for each sample
#' @param method (character) One of "Beauchamp" or "Myers" indicating the method
#'                           "Beauchamp": Average the expression of multiple
#'                                        probes per gene
#'                           "Myers": Three step probe selection based on 
#'                                    mean intensity and probe-probe correlation
#' @param verbose (logical) Verbosity option
#'
#' @return (list) A list containing three elements:
#'                Donor (character) 
#'                    Donor name
#'                GeneExpression (tibble)
#'                    Tibble containing gene-by-sample expression values
#'                SampleInfo (tibble) 
#'                    Tibble containing metadata for each sample
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
    
    #Select probe per gene based on maximal average correlation to all other
    #probes for that gene
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