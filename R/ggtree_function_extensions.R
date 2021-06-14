get_clade_nodes <- function (x, y, z){ # give phylogroups, clusters,tree


  ### cycle through the clusters - get the mrca, then subtree out to put this information on the tree:
  ###make df for the nodes, with the labels as well
  cluster_nodes <- data.frame(matrix(ncol = 3, nrow = 0))
  colnames(cluster_nodes) <- c("node","poppunk_cluster","phylogroup")

  need_to_break = F


  for (c in unique(y$Cluster)) {
    temp_df <- subset.data.frame(y, Cluster == c)

    ## get the phylogroup that each poppunk cluster belongs to:
    strains_in_poppunk_cluster <- data.frame(rownames(temp_df))
    colnames(strains_in_poppunk_cluster) <- c("strain") # rename columns
    ## now merge in to the strains_and_phylgroups_df to get the phylogroup, then check to make sure there is only one phylogroup:
    check_df <- merge(x, strains_in_poppunk_cluster, by.x = "strain", by.y = "strain")
    #check number of phylogroups in this:
    if (length(unique(check_df$phylogroup)) > 1) {
      print(paste0("Error - more than one phylogroup in cluster ",c,"!", sep = ""))
      need_to_break = TRUE
    }
    ##check if NA:
    if (is.na(unique(check_df$phylogroup)[1]) == TRUE) {
      print(paste0("Error - phylogroup in cluster ",c," is NA!", sep = ""))
      need_to_break = TRUE
    }
    ### now break here - if needed
    if (need_to_break == TRUE) {
      break
    }

    ## if here, then only one phylogroup in the list:
    this_phylogroup <- unique(check_df$phylogroup[1])
    # add this to the dataframe in the loop below


    # get the mrca to annotate the tree with the cluster
    ##if (nrow(temp_df) > 2) { # use if want poppunk clusters with more than 2:

    mrca_node <- ape::getMRCA(z, rownames(temp_df))
    # if only one assembly in the poppunk cluster, then the MRCA will be NULL - so need to get this node:
    if (is.null(mrca_node)) {
      # get node and add to the to_rbind dataframe
      mrca_node <- subset.data.frame(tibble::as_tibble(z), label == check_df$strain)$node
    }

    to_rbind <- data.frame(mrca_node, c,this_phylogroup)
    colnames(to_rbind) <- colnames(cluster_nodes)
    cluster_nodes <- rbind(cluster_nodes, to_rbind)

    #}
  }

  # ### now continue here - if needed
  # if (need_to_break == TRUE) {
  #   continue
  # }

  clade_labels <- cluster_nodes$node
  names(clade_labels) <- cluster_nodes$poppunk_cluster


  return(clade_labels)
}


add_clades_to_tree <- function(drawn_tree, clade_labels,highlight_clades,aligned,extend_to_value,highlight_alpha,label_offset) {

  #1 - tree already drawn as ggtree object
  #2 - list of the internal nodes to be labelled / highlighted as clades , with the desired labels as names of the list / vector
  #3 - TRUE/FALSE as to whether to include the highlight blocks as well as the line cladelabel

  # drawn_tree <- t2
  # clade_labels <- clade_labels
  # highlight_clades <- TRUE
  # aligned <- T


  temp_ggplot <- drawn_tree

  ## loop though and add the clade labels:
  for (l in 1:length(clade_labels)) {
    temp_ggplot <- temp_ggplot +
      ggtree::geom_cladelabel(node=clade_labels[l], label=names(clade_labels)[l], offset = label_offset, align = aligned)
  }

  is.odd <- function(x) x %% 2 != 0

  #arrange the clade highlighting from top of plot to lower , i.e. from high to low for y values for the nodes in clade_labels:
  col_df <- temp_ggplot$data %>%
    right_join(data.frame(node = clade_labels)) %>%
    arrange(-y) %>%
    mutate(color = unlist(purrr::map(seq_along(y), ~ ifelse(is.odd(.x),"#a9a9a9","#d9d9d9"))))


  if (isTRUE(highlight_clades)) {

    for (i in 1:nrow(col_df)) {

      temp_ggplot <- temp_ggplot +
        ggtree::geom_hilight(node=col_df[i,]$node, fill = col_df[i,]$color, extendto = max(temp_ggplot$data$x) + extend_to_value + 0.001 , alpha = highlight_alpha)
    }
  }

  ## try geom_hilight:
  # for (l in 1:length(clade_labels)) {
  #   ## check if odd:
  #   if(is.odd(l)) {
  #     col = "#a9a9a9"   extend = max(col_df$x) - col_df[i,]$x
  #   } else{
  #     col = "#d9d9d9"
  #   }
  #   ## now add to plot
  #   temp_ggplot <- temp_ggplot +
  #     geom_hilight(node=clade_labels[l], fill = col, extend = 0.2, alpha = 0.25)
  #
  # }




  return(temp_ggplot)


}


remove_tips_from_tree <- function(tree,tips) {

  require(ggtree)

  temp_tree <- tree

  for (tipname in tips) {
    temp_tree <- drop.tip(temp_tree, tipname)
  }

  return(temp_tree)

}


collapse_tree <- function(tree,nodes,type,colours) { # provide list of nodes

  temp_tree <- tree

  if ( ! (type  %in% c("min","max","mixed"))) {
    stop("Please provide either 'min', 'max', or 'mixed'")
  }

  if ( isFALSE(colours)) {

    for (l in 1:length(nodes)) {
      node_num = nodes[l]
      temp_tree <- ggtree::collapse(temp_tree, node = node_num, type)
    }
  } else {

    for (l in 1:length(nodes)) {
      node_num = nodes[l]
      hex_colour <- colours[l]
      temp_tree <- ggtree::collapse(temp_tree, node = node_num, type, fill= hex_colour, color = hex_colour, alpha = .5)
    }


  }


  return(temp_tree)
}







