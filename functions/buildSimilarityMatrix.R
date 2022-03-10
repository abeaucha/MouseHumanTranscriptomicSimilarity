#' Build a similarity matrix
#'
#' @description
#' `buildSimilarityMatrix` computes the similarity between inputs x1 and x2
#'
#' @details
#' 
#'
#' @param x1 A matrix or data tree containing values to use to calculate similarity.
#' @param x2 A matrix or data tree containing values to use to calculate similarity.
#' @param method A character scalar indicating which similarity metric to use. 
#'               One of "correlation", "cosine", "euclidean".
#' @param treefield A character scalar indicating which tree field contains the values to use.
#'                  Only used when inputs are data trees.
buildSimilarityMatrix <- function(x1, x2, method = "correlation", treefield = NULL){
  
  
  #Check if x1 is either a matrix or data tree
  if (!isTree(x1) & !is.matrix(x1)) {
    stop(paste("Argument `x1` must be a matrix or data tree. Got class", paste(class(x1), collapse = "; ")))
  }
  
  #Check if x2 is either a matrix or data tree
  if (!isTree(x2) & !is.matrix(x2)) {
    stop(paste("Argument `x2` must be a matrix or data tree. Got class", paste(class(x2), collapse = "; ")))
  }
  
  #Condition for data tree inputs
  if (isTree(x1) & isTree(x2)) {
    
    #Check for treefield existence
    if (is.null(treefield)) {
      stop("No value passed to argument `treefield`. Necessary for data trees.")
    } else {
      if (!(treefield %in% x1[["fields"]])){
        stop(paste("`treefield` value", treefield, "not found in `x1`"))
      }
      if (!(treefield %in% x2[["fields"]])){
        stop(paste("`treefield` value", treefield, "not found in `x2`"))
      }
    }
    
    #Build matrices from leaf nodes in tree
    mat1 <- x1$Get(treefield, filterFun = data.tree::isLeaf)
    mat2 <- x2$Get(treefield, filterFun = data.tree::isLeaf)
    
    #Center data clouds in gene space
    #Probably too specific for this function
    mat1 <- mat1 - rowMeans(mat1)
    mat2 <- mat2 - rowMeans(mat2)
    
  } else if (is.matrix(x1) & is.matrix(x2)) {
    
    mat1 <- x1
    mat2 <- x2
    
  } else {
    stop("Arguments `x1` and `x2` must both be either matrices or data trees")
  }
  
  #Bind matrices into one
  mat <- cbind(mat1, mat2)
  
  #Indices to extract similarity matrices
  nCol1 <- ncol(mat1)
  ind1 <- 1:nCol1
  ind2 <- (nCol1+1):ncol(mat)
  
  #Define similarity metric
  if (method == "correlation"){
    similarity <- cor
  } else if (method == "cosine") {
    similarity <- cosineSimilarity
  } else if (method == "euclidean"){
    similarity <- distEuclidean
  } else {
    stop(paste("Invalid value passed to argument `method`:", method))
  }
  
  #Compute similarity
  matSimAll <- similarity(mat)
  
  #Extract cross-species similarity
  matSim <- matSimAll[ind1, ind2]
  
  return(matSim)
}