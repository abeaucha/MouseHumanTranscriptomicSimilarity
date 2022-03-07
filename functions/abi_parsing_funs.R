parse_abi_hierarchy <-
  function(abi){
    abi_list <- fromJSON(file = abi)$msg[[1]]
    tree <- recursive_build_tree(abi_list)
  }

write_abi_hierarchy <-
  function(tree, file = NULL){
    process_node <-
      function(node){
        node_list <- c(node$name, lapply(node$fields, function(a) getElement(node, a)))
        names(node_list) <- c("name", node$fields)
        node_list$children <- lapply(node$children, process_node)
        
        node_list
      }
    
    tree_list <- list(msg = list())
    tree_list$msg <- list(process_node(tree))

    json <- toJSON(tree_list)

    if(!is.null(file))
      cat(json, file = file)

    return(json)
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

filter_keep_parents <-
  function(tree, fun, keep_children = FALSE){
    node <- Clone(tree)
    
    prune_children <- function(node, fun, parent){
      if(is.null(node$children)){
        if(!fun(node)) parent$RemoveChild(node$name)
      } else {
        node_matches <- fun(node)
        
        if(keep_children && node_matches){
          return(NULL)
        } else {
          lapply(node$children, prune_children, fun, node)
          if(!node_matches && length(node$children) == 0){
            parent$RemoveChild(node$name)
          }
        }
      }
    }
    
    lapply(node$children, prune_children, fun, node)
    
    if(is.null(node$children)){
      return(NULL)
    } else {
      return(node)
    }
  }

add_structure <-
  function(tree, name, parent = NULL, id = NULL, filterFun = NULL, lateralized = TRUE, colour = NULL){
    parent_passed <- is.null(parent)
    filter_passed <- is.null(filterFun)
    
    if(parent_passed == filter_passed)
      stop("Only one of parent and filter fun must be passed (not both)")
    
    parent_list <- Traverse(tree, filterFun = function(n) n$name == parent)
    if(length(parent_list) != 1) 
      stop("The parent node isn't uniquely specified, you may need to use a filterFun")
    
    if(is.null(id))
      id <- get_new_id(tree)
    
    parent_node <- parent_list[[1]]
    parent_path <- parent_node$structure_id_path
    parent_id <- parent_node$id
    parent_colour <- parent_node$color_hex_triplet
    
    if(is.null(colour))
      colour <- nudge_colour(parent_colour)
    
    new_node <-
      Node$new(name = name, id = id, structure_id_path = paste0(parent_path, "/", id), parent_id = parent_id
               , color_hex_triplet = colour)
    
    parent_node$AddChildNode(new_node)
    
    if(lateralized){
      add_structure(tree, paste0(name, ", left"), name, lateralized = FALSE, colour = colour)
      add_structure(tree, paste0(name, ", right"), name, lateralized = FALSE, colour = colour)
    }
    
    invisible(tree)
  }

get_new_id <-
  function(tree, chunk_size = 100000){
    labels <- tree$Get("id")
    min_label <- min(labels)
    
    
    potential_labels <- min_label:(min_label + chunk_size)
    usable <- setdiff(potential_labels, labels)
    
    min(usable)
  }

nudge_colour <-
  function(colour, percent = .1){
    c(substr(colour, 1,2)
      , substr(colour, 3,4)
      , substr(colour, 5,6)) %>%
      paste0("0x", .) %>%
      strtoi %>%
      `+`(255 * (runif(3) * (2 * percent) - percent)) %>% ##Add or subtract between -percent and percent% of the color val
      round %>%
      clamp(0,255) %>%
      as.hexmode %>%
      paste0(collapse = "") %>%
      toupper
  }

clamp <- 
  function(num, min, max)
    sapply(num, function(n){
      if(n < min) return(min)
      if(n > max) return(max)
      return(n)
    })



human_volumes_to_hierarchy <-
  function(tree, pond_volume_summary, defs = "AAL_allen_coarse_merge_after_edits.csv"){
    tree <- Clone(tree)
    
    labels <- read.csv(defs)
    unilats <- labels %>% filter(Allen_left == Allen_right)
    bilats <- labels %>% filter(Allen_left != Allen_right)
    
    bilats_long <- 
      bilats %>% 
      gather(hemisphere, label, Allen_left, Allen_right) %>%
      mutate(hemisphere = 
               sub("Allen_", "", hemisphere) %>% 
               sub("^r", "R", .) %>% 
               sub("^l", "L", .)
             , Structure = paste0(hemisphere, "_", RMINC:::fix_names(Structure)))
    
    unilats_long <-
      unilats %>%
      rename(label = Allen_left) %>%
      select(-Allen_right) 
    
    labels_long <- bind_rows(bilats_long, unilats_long)
    
    gyral_volumes <- 
      select(pond_volume_summary, scan, matches("lobeVolume")) %>%
      mutate(row = row_number()) %>%
      gather(structure, volume, -scan, -row) %>%
      mutate(structure = sub("lobeVolume_", "", structure)) %>%
      spread(structure, volume)
    
    subcortical_volumes <-
      select(pond_volume_summary, scan, Vermal_I_II:right.ventral.posterior.nucleus) %>%
      mutate(row = row_number()) %>%
      gather(structure, volume, -scan, -row) %>%
      mutate(structure = 
               gsub("\\.+", "_", structure) %>% 
               gsub("^_|_$", "", .) %>%
               sub("right", "Right", .)  %>% 
               sub("left", "Left", .)) %>%
      spread(structure, volume)
    
    volume <-
      inner_join(gyral_volumes, subcortical_volumes, by = c("scan", "row"))
    
    tree$Do(function(n){
      id <- n$id
      if(id %in% labels_long$label){
        if(!isLeaf(n)) stop("Something went wrong, trying to assign a value to a non-leaf")
        
        structure <- labels_long$Structure[labels_long$label == id]
        print(structure)
        n$volumes <- volume[,structure]
        n$meanVolume <- mean(n$volumes,na.rm = TRUE)
      }
    })
    
    tree$Do(function(n) n$meanVolume <- Aggregate(n, "meanVolume", sum), traversal = "post-order")
    tree$Do(function(n) n$volumes <- Aggregate(n, "volumes", rowSums)
            , traversal = "post-order"
            , filterFun = isNotLeaf)
      
    tree
  }


