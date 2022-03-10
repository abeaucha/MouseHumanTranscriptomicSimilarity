#' Compute the cosine similarity
#'
#' @param x1 (matrix) A numeric matrix 
#' @param x2 (matrix) A numeric matrix
#'
#' @return A square matrix containing the cosine similarity between the columns of m1 (and m2)
cosineSimilarity <- function(x1, x2=NULL){
  
  if(is.null(x2)){
    x2 <- x1
  }
  
  numerator <- t(x1) %*% x2
  denominator1 <- as.matrix(diag(t(x1) %*% x1))
  denominator2 <- as.matrix(diag(t(x2) %*% x2))
  denominator <- sqrt(denominator1 %*% t(denominator2))
  
  return(numerator/denominator)
  
}

#' Compute the Euclidean distance
#'
#' @param x A numeric matrix
#'
#' @return A square matrix containing the Euclidean distance between all columns
distEuclidean <- function(x){as.matrix(dist(t(x), method = "euclidean"))}

#' Extract a sub-matrix
#'
#' @param x 
#' @param rowInd 
#' @param colInd 
#'
#' @return
extractSubMat <- function(x, rowInd, colInd = NULL){
  if(is.null(colInd)){colInd <- rowInd}
  return(x[rowInd, colInd])
}