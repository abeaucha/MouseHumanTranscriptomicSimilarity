#' Intersect gene homologs
#'
#' @param data 
#' @param homologs 
#'
#' @return
intersectGeneHomologs <- function(data, homologs){
  
  require(data.tree)
  require(assertive)
  
  #Check that data is a named list
  assert_is_list(data)
  assert_has_names(data)
  
  #Check that names are Mouse and Human
  if(!(any(c("Mouse", "Human") %in% names(data)))){
    stop("`data` should have elements named `Mouse` and `Human`",
         call. = FALSE)
  }
  
  #If homologs is a character, load the file
  if(is.character(homologs)){
    homologs <- suppressMessages(read.csv(homologs, 
                                          stringsAsFactors = FALSE)) 
  }
  
  #Check that homologs is a data frame
  assert_is_data.frame(homologs)
  
  #Check that homologs has columns Mouse and Human
  if(!(any(c("Mouse", "Human") %in% colnames(homologs)))){
    stop("`homologs` should have columns `Mouse` and `Human`",
         call. = FALSE)
  }
  
  #Subset those columns
  homologs <- homologs[, c("Mouse", "Human")]
  
  #Extract mouse and human data
  mouse <- data[["Mouse"]]
  human <- data[["Human"]]
  
  #Identifying homologs in the mouse and human data
  indMouse <- homologs$Mouse %in% mouse$Gene
  indHuman <- homologs$Human %in% human$Gene
  
  #Identify the intersection of all data frames
  homologs <- homologs[indMouse & indHuman, ]
  
  
  #Control flow for data frames and trees
  if(is.data.frame(mouse) & is.data.frame(human)){
    
    #Subset the mouse and human data to be in the intersection
    mouse <- mouse[mouse$Gene %in% homologs$Mouse, ]
    human <- human[human$Gene %in% homologs$Human, ]
    
    #Rename the mouse genes using human symbols
    mouse$Gene <- homologs[match(mouse$Gene, homologs$Mouse), ]$Human
    
    #Order the data frames by gene
    mouse <- mouse[order(mouse$Gene), ]
    human <- human[order(human$Gene), ]
    
  } else if(isTree(mouse) & isTree(human)){
    
    mouse_new <- Clone(mouse) 
    human_new <- Clone(human)
    
    mouse_new$Do(function(node){
      
      #Subset data to be in the intersection
      indHom <- node$Gene %in% homologs$Mouse
      node$Gene <- node$Gene[indHom]
      node$Expression <- node$Expression[indHom]
      
      #Rename genes
      node$GeneMouse <- node$Gene
      node$Gene <- homologs[match(node$Gene, homologs$Mouse), ]$Human
      
      #Order by human names
      indOrder <- order(node$Gene)
      node$Gene <- node$Gene[indOrder]
      node$GeneMouse <- node$GeneMouse[indOrder]
      node$Expression <- node$Expression[indOrder]
    })
    
    human_new$Do(function(node){
      indHom <- node$Gene %in% homologs$Human
      node$Gene <- node$Gene[indHom]
      node$Expression <- node$Expression[indHom]
      
      node$GeneMouse <- homologs[match(node$Gene, homologs$Human), ]$Mouse
      
      indOrder <- order(node$Gene)
      node$Gene <- node$Gene[indOrder]
      node$GeneMouse <- node$GeneMouse[indOrder]
      node$Expression <- node$Expression[indOrder]
    })
    
    mouse <- Clone(mouse_new)
    human <- Clone(human_new)
    
  } else {
    stop("`data` contains neither data frames nor data trees.",
         call. = FALSE)
  }
  
  #Return
  out <- list(Mouse = mouse,
              Human = human)
  
  return(out)
  
}


isTree <- function(x){
  any(class(x) == "Node") & 
    any(class(x) == "R6")
}



imputer <- function(data, strategy = "mean", axis = "columns"){
  
  if (strategy == "mean"){
    impute <- mean
  } else if (strategy == "median"){
    impute <- median
  }
  
  if (axis == "rows"){
    
    for(i in 1:nrow(data)){
      if(any(is.na(data[i,]))){
        data[i,is.na(data[i,])] <- impute(data[i,], na.rm = TRUE)
        }
    }
    
  } else {
    
    for (i in 1:ncol(data)){
      if(any(is.na(data[,i]))){
        data[is.na(data[i,]),i] <- impute(data[,i], na.rm = TRUE)
        }
    }
    
  }
  
  return(data)
}


scaler <- function(data, scale = TRUE, axis = "columns"){
  
  if (axis == "columns"){
    
    colMeansMat <- matrix(colMeans(data), 
                          nrow = nrow(data), 
                          ncol = ncol(data), 
                          byrow = T)
    out <- data - colMeansMat
    
    if (scale == TRUE){
      
      sigma <- sqrt(colSums(out^2)/(nrow(data)-1))
      sigmaMat <- matrix(sigma, 
                         nrow = nrow(data), 
                         ncol = ncol(data), 
                         byrow = T)
      out <- out/sigmaMat
      
    } 
    
  } else if (axis == "rows") {
    
    rowMeansMat <- matrix(rowMeans(data), 
                          nrow = nrow(data), 
                          ncol = ncol(data), 
                          byrow = F)
    out <- data - rowMeansMat
    
    if (scale == TRUE){
      
      sigma <- sqrt(rowSums(out^2)/(ncol(data)-1))
      sigmaMat <- matrix(sigma, 
                         nrow = nrow(data), 
                         ncol = ncol(data), 
                         byrow = F)
      out <- out/sigmaMat
      
    }
  }
  
  return(out)
  
} 