#' Check if object is a data tree
#'
#' @description 
#' `isTree` checks whether its argument is a data tree
#'
#' @param x Object to test for tree structure
#' 
#' @return Boolean value TRUE or FALSE
isTree <- function(x){any(class(x) == "Node") & any(class(x) == "R6")}


parse_abi_hierarchy <-
  function(abi){
    abi_list <- fromJSON(file = abi)$msg[[1]]
    tree <- recursive_build_tree(abi_list)
  }

recursive_build_tree <- function(tree_list){
  children_list <- tree_list$children
  node_list <- tree_list[names(tree_list) != "children"]
  
  node <- do.call(Node$new, node_list)
  
  if(is.null(children_list) || length(children_list) == 0){
    return(node)
  } else {
    lapply(children_list, function(cl){
      child_node <- recursive_build_tree(cl)
      node$AddChildNode(child_node)
      
      NULL
    })
    
    return(node)
  }
}

#' Attach cut points to a tree node
#' 
#' @description 
#' This function attaches cut points to a data tree node. The cut
#' points are used by the pruning functions pruneAtNode and 
#' pruneBelowNode to prune the tree based on those nodes.
#' 
#' @param node A data tree node 
#' @param cutPoints A character vector containing the names of tree
#'                  nodes to use as cut points
#'                  
#' @return NULL                
attachCutPoints <- function(node, cutPoints){node[["cutPoints"]] <- cutPoints}



#' Prune data tree at node
#' 
#' @description 
#' This pruning function evaluates whether the name of the node is
#' in the cut points (stored in the node field "cutPoints"). When 
#' passed to the `pruneFun` argument in the `Prune` function from 
#' `data.tree` it will prune the data tree at the nodes specified 
#' in the "cutPoints" field.
#'
#' @param node A data tree node
#' 
#' @return Boolean value TRUE or FALSE 
pruneAtNode <- function(node){
  
  if(is.null(node[["cutPoints"]])){
    stop(paste("Node", node[["name"]], "has no field 'cutPoints'. Cut points can be attached using the function attachCutPoints()."))
  }
  
  return(!(node$name %in% node$cutPoints))
}



#' Prune data tree below node
#' 
#' @description 
#' This pruning function evaluates whether the node is a child of 
#' the cut points (stored in the node field "cutPoints"). When 
#' passed to the `pruneFun` argument in the `Prune` function from 
#' `data.tree` it will prune the data tree below the nodes specified
#' in the "cutPoints" field.
#'
#' @param node A data tree node
#' 
#' @return Boolean value TRUE or FALSE 
pruneBelowNode <- function(node){
  
  if(is.null(node[["cutPoints"]])){
    stop(paste("Node", node[["name"]], "has no field 'cutPoints'. Cut points can be attached using the function attachCutPoints()."))
  }
  
  return(!any(node$cutPoints %in% node$path[-length(node$path)]))
  }



#' Prune a tree based on node names
#' 
#' @description 
#' 
#' @param anatTree The data tree to prune
#' @param nodes A character vector of node names to use to prune the
#'              tree
#' @param method One of "AtNode" or "BelowNode" indicating whether
#'               to prune the tree at or below the nodes specified              
#' 
#' @return 
pruneAnatTree <- function(anatTree, nodes, method = "AtNode"){
  
  require(data.tree)
  
  #Select where to prune
  if(method == "AtNode"){
    pruningFunction <- pruneAtNode
  } else if (method == "BelowNode"){
    pruningFunction <- pruneBelowNode
  } else {
    stop(paste("Argument `method` must be one of 'AtNode' or 'BelowNode'. Got", method))
  }
  
  #Attach cut points to tree
  anatTree$Do(attachCutPoints, cutPoints = nodes)
  
  #Prune cut points from tree 
  Prune(anatTree, pruningFunction)
  
  outMessage <- paste("Tree pruned",
                      ifelse(method == "AtNode", "at nodes", "below nodes"),
                      paste(nodes, collapse = ", "))
  
  return(outMessage) 
}



## Function: hanatToAtlas -------------------------------------------
hanatToAtlas <- function(anatTree, labelVolume){
  
  require(RMINC)
  require(data.tree)
  
  anatTree$Do(function(node){
    node$label.new <- min(node$label)
  }, traversal = "post-order")
  
  out <- hanatToVolume(anatTree, labelVolume, "label.new")
  out <- as.numeric(out)
  class(out) <- class(labelVolume)
  attributes(out) <- attributes(labelVolume)
  attr(out, "dim") <- NULL
  
  return(out)
}


## Function: hanatToAtlasDefs ------------------------------------
#
# Generates new atlas definitions based on a pruned anatomical 
# tree.
#
# Arguments: 
# anatTree: The pruned tree for which to generate an atlas. Must
#           have fields `label` and `label.new`.
#
# Return:
# A data frame giving the new label definitions
# -------------------------------------------------------------------
hanatToAtlasDefs <- function(anatTree){
  
  require(data.tree)
  
  if(is.null(anatTree[["label.new"]])){
    anatTree$Do(function(node){
      node$label.new <- min(node$label)
    }, traversal = "post-order")
  }
  
  labelsNew <- anatTree$Get("label.new", filterFun = isLeaf)
  
  return(data.frame(Structure = names(labelsNew),
                    Label = labelsNew,
                    row.names = NULL,
                    stringsAsFactors = FALSE))
}