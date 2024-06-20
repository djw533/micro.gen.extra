#' Find internal clade nodes ancestral to named/numbered sets of tip labels within a second layer of clustering (within a tree)
#'
#' For a tree object with an associated set of clustering (e.g. phylogroups defined by fastANI), and a finer set
#' of clustering (e.g. multi-level clustering of bacterial genomes using fastbaps/PopPUNK), where all of the clusters
#' in the latter clustering are sets entirely within the former clusters, retrieve the mrca for those sets of tips
#' alongside the fastbaps/PopPUNK cluster number, as a named vector list.
#'
#' @param x Dataframe with two columns: "strain","phylogroup", where strain corresponds with tip-label
#' @param y Dataframe with two columns: "strain", "Cluster", where strain corresponds with tip-label
#' @param z Tree object, where the tip-labels correspond with "strain" in x and y
#'
#' @return Named vector list where the values are the internal node numbers, and the names are the cluster numbers.
#'
#' @examples
#' get_within_species_cluster_mrcas(x = phylogroups_dataframe, y = clusters_dataframe, z = ggtree_object)

get_within_species_cluster_mrcas <- function (x, y, z){ # give phylogroups, clusters,tree

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



#' Find internal clade nodes ancestral to named/numbered sets of tip labels (within a tree)
#'
#' For a tree object with an associated set of clustering (e.g. phylogroups defined by fastANI, or clustering using
#' fastbabs/PopPUNK), retrieve the mrca for those sets of tips alongside the fastbaps/PopPUNK cluster number,
#' as a named vector list.
#'
#' @param tree Tree object, e.g. read in using ape::read.tree()
#' @param clusters_df Dataframe with two columns: "strain", "cluster", where strain corresponds with tip-label
#'
#' @return Named vector list where the values are the internal node numbers, and the names are the cluster numbers.
#'
#' @examples
#' get_clade_nodes(tree_object, clusters_dataframe)

get_clade_nodes <- function(tree, clusters_df) {


  #check the colnames are right:

  if (isFALSE("strain" %in% colnames(clusters_df)) | isFALSE("cluster" %in% colnames(clusters_df))) {
    stop(glue::glue('Please make sure that "strain" and "cluster" are colnames in clusters_df'))
  }



  #ggtree object to extract out nodes for clusters that are singletons :
  tree_df <- ggtree::ggtree(tree)$data


  #get the mrca nodes
  clusters_vector <- clusters_df %>%
    dplyr::group_by(cluster) %>%
    dplyr::mutate(strain_concatenated = list(strain)) %>% # get lists of the tips for each cluster
    dplyr::select(-strain) %>%
    unique() %>% # remove duplicates after getting rid of the strain
    #dplyr::mutate(mrca_node = purrr::map(strain_concatenated, ~ ape::getMRCA(tree, .x))) %>% # get the mrca node from the list of tips in each cluster
    dplyr::mutate(mrca_node = purrr::map(strain_concatenated, ~ ifelse(length(.x) > 1, ape::getMRCA(tree, .x), tree_df %>% filter(label == .x, isTip == T) %>% pull(node) ))) %>% # get the mrca node from the list of tips in each cluster
    pull(mrca_node, name = cluster) #pull out named vector


  return(clusters_vector)

}



#' Highlight clades over a ggtree object given a named vector of internal nodes and labels to draw
#'
#' Given a named vector of internal tree nodes and labels (e.g. cluster number / name), move through each cluster name
#' and apply geom_hilight and geom_cladelabel from ggtree onto the tree, with alternating shades of grey from top
#' to bottom.
#'
#' @param drawn_tree ggtree object for labels to be added to
#' @param clade_labels Named vector of internal node and label (as names), e.g. output from micro.gen.extra::get_clade_nodes()
#' @param highlight_clades TRUE/FALSE as to whether to include the ggtree::geom_hilight [Default = TRUE]
#' @param aligned TRUE/FALSE as to whether the ggtree::geom_hilights should align [Default = TRUE]
#' @param extend_to_value Value to extend the ggtree::geom_hilight to. Float. [Default = 0.0]
#' @param highlight_alpha Value to set the alpha of the ggtree::geom_hilight. Set between 0 and 1. [Default = 0.4]
#' @param label_offset Value to offset the labels of the clade_labels. [Default = F]
#' @param gradient T/F for whether hilighted box should have a gradient (one of "rt" or "tr" [root-to-tip or tip-to-root]) [Default = "tr"]
#' @param gradient_drcton Gradient direction [Default = 0.0]
#'
#' @return ggtree object with clades (corresponding to the internal nodes supplied) have been highlighted and labelled.
#'
#' @examples
#' add_clades_to_tree(ggtree_object, clade_labels_vector)


add_clades_to_tree <- function(drawn_tree,
                               clade_labels,
                               highlight_clades = TRUE,
                               aligned = TRUE,
                               extend_to_value = 0.0,
                               highlight_alpha = 0.4,
                               label_offset = 0.0,
                               gradient = F,
                               gradient_direction = "tr") {

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
    right_join(data.frame(node = unlist(clade_labels))) %>%
    arrange(-y) %>%
    mutate(color = unlist(purrr::map(seq_along(y), ~ ifelse(is.odd(.x),"#a9a9a9","#d9d9d9"))))


  if (isTRUE(highlight_clades)) {

    for (i in 1:nrow(col_df)) {

      if (isFALSE(gradient)) {
        temp_ggplot <- temp_ggplot +
          ggtree::geom_hilight(node=col_df[i,]$node,
                               fill = col_df[i,]$color,
                               extendto = max(temp_ggplot$data$x) + extend_to_value + 0.001 ,
                               alpha = highlight_alpha)
      } else if (isTRUE(gradient)) {
        temp_ggplot <- temp_ggplot +
          ggtree::geom_hilight(node=col_df[i,]$node,
                               fill = col_df[i,]$color,
                               extendto = max(temp_ggplot$data$x) + extend_to_value + 0.001 ,
                               alpha = highlight_alpha,
                               type = "gradient",
                               gradient.direction = gradient_direction)
      }
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



#' Prune a list of tips from a tree
#'
#' Given a list of tips, remove these from a tree object and return the pruned tree
#'
#' @param tree Tree object
#' @param tips List of tips as present in tree
#'
#' @return Tree object with requested tips pruned out
#'
#' @examples
#' remove_tips_from_tree(tree_object, list_of_tips)

remove_tips_from_tree <- function(tree,tips) {

  require(ape)

  temp_tree <- tree

  for (tipname in tips) {
    temp_tree <- ape::drop.tip(temp_tree, tipname)
  }

  return(temp_tree)

}


#' Collapse a set of nodes in a ggtree object
#'
#' Given a list of nodes, collapse these nodes in a ggtree object
#'
#' @param tree Tree object
#' @param nodes List of internal nodes to collapse
#' @param type Type of collapsing to run - "min", "max" or "mixed". [Default = "max"]
#' @param colours List of hex colours to colour the collapsed cones, order of the colours corresponds
#' with the list of nodes given [Default = "none"]
#'
#' @return Tree object with requested tips pruned out
#'
#' @examples
#' remove_tips_from_tree(tree_object, list_of_tips)


collapse_tree <- function(tree,nodes,type = "max",colours = "none") {

  temp_tree <- tree

  if ( ! (type  %in% c("min","max","mixed"))) {
    stop("Please provide either 'min', 'max', or 'mixed'")
  }

  if ( colours == "none") {

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



#' Put a multi-column gheatmap onto a ggtree object with gaps between columns
#'
#' Put a multi-column gheatmap onto a ggtree object with gaps between columns
#'
#' @param ggtree_object Tree object
#' @param gheatmap_df Dataframe to be put along tree as gheatmap. Rownames should correspond with the tip labels in the tree.
#' @param col_width Column width for each column in the heatmap [Default = 1/5 width of the tree]
#' @param offset_num Initial offset from the tree [Default = 0]
#' @param gap Gap between each column [Default = 1/10 of value of column width]
#'
#' @return Tree object with gheatmap added, with the gaps
#'
#' @examples
#' remove_tips_from_tree(tree_object, list_of_tips)


gheatmap_with_gaps <- function(ggtree_object,
                               gheatmap_df,
                               col_width = max(ggtree_object$data$x)/5,
                               offset_num = 0,
                               gap = col_width / 10) {

  iteration <- 0
  temp_plot <- ggtree_object

  #check the rownames are all in the tree:

  if ( isFALSE(all(rownames(gheatmap_df) %in% subset(ggtree_object$data, isTip == TRUE)$label)) ) {
    stop("All rownames in the gheatmap should be tips in the ggtree object. Stopping")
  }


  for (i in colnames(gheatmap_df)) {

    if (iteration == 0) {
      offset_num <- 0
    }

    iteration = iteration + 1

    temp_data <- gheatmap_df %>%
      select(i)

    temp_plot <- temp_plot %>% ggtree::gheatmap(temp_data,
                                                colnames_position = "top", colnames_offset_y = 15,
                                                color = NULL,
                                                width = col_width,
                                                offset = offset_num) +
      #scale_fill_manual(values = c(fill_colours)) +
      ylim(0,max(temp_plot$data$y) + 20)


    offset_num <- (iteration * (max(temp_plot$data$x) * col_width)) + (gap * iteration)

  }

  return(temp_plot)
}


t6_gheatmap_with_gaps <- function(ggtree_object,
                                  t6_subtypes_df,
                                  col_width = max(ggtree_object$data$x)/5,
                                  colname_offset = 0,
                                  init_col_offset = 0,
                                  col_gap = col_width / 10) {


  cols <- c("predicted_T6SSs" = "blue",
            "i1" = "#3cb44b",
            "i2" = "#ffe119",
            "i3" = "#e6194b",
            "i4a" = "#4363d8",
            "i4b" = "#ff1493",
            "i5" = "#911eb4",
            "i3v3" = "#76b7b2",
            "i3v1" = "#ff9da7",
            "i3v2" = "#f28e2b",
            "Undetermined" = "grey",
            "rejected_T6SSs" = "grey",
            "contig_broken_T6SS" = "grey")


  iteration <- 0
  temp_plot <- t6_subtypes_df %>%
    tibble::rownames_to_column(var = "strain")

  #check the rownames are all in the tree:

  if ( isFALSE(all(rownames(gheatmap_df) %in% subset(ggtree_object$data, isTip == TRUE)$label)) ) {
    stop("All rownames in the gheatmap should be tips in the ggtree object. Stopping")
  }

  for (i in colnames(gheatmap_df)) {

    if (iteration == 0) {
      offset_num <- init_col_offset
    }

    iteration = iteration + 1

    temp_data <- gheatmap_df %>%
      select(i)

    temp_plot <- temp_plot %>% ggtree::gheatmap(temp_data,
                                                colnames_position = "top",
                                                colnames_offset_y = 15,
                                                colnames_angle = 90,
                                                hjust = 0,
                                                color = NULL,
                                                width = col_width,
                                                offset = colname_offset) +
      scale_fill_gradient(low = "white", high = cols[i]) +
      ggnewscale::new_scale_fill() +
      ylim(0,max(temp_plot$data$y) + 20)


    offset_num <- (iteration * (max(temp_plot$data$x) * col_width)) + (col_gap * iteration)

  }


  return(temp_plot)


}


#' Subsample a phylogenetic tree
#'
#' Subsample a phylogenetic tree by randomly selecting one node for all clades/subtrees of a maximum tip-root distance
#'
#' @param input_tree Tree object (e.g. from ape)
#' @param cut_dist max tip-root distance to select clades to be sub-sampled
#' @param plot_clades_to_prune plot original input tree with clades to be pruned annotated [Default = False]
#'
#' @return Pruned tree object
#'
#' @examples
#' subsample_tree(input_tree = tree_object, cut_dist = 0.02)

subsample_tree <- function(input_tree, cut_dist, plot_clades_to_prune = F) {

  #2 take all subtrees in tree
  tree.subtrees <- ape::subtrees(input_tree)

  #3 filter to the distance
  subtree_distances <- data.frame(subtree_number = seq_along(tree.subtrees)) %>%
    mutate(tips = purrr::map(subtree_number, ~ sort(tree.subtrees[[.x]]$tip.label))) %>%
    mutate(max.dist = as.numeric(purrr::map(subtree_number, ~ max(ape::cophenetic.phylo(tree.subtrees[[.x]])))))

  #4 get all nested subtrees and remove
  subtree_distances.filtered <- subtree_distances %>%
    filter(max.dist <= cut_dist) %>%
    #check if they are nested
    mutate(nested = unlist(
      purrr::map2(tips, subtree_number, function(a,x)
        any(unlist(    # does subtree b contain all of the tips that are in any other subtree in the dataset?
          purrr::map2(tips, subtree_number, function(b,y)
            all(a %in% b) & x != y) # return whether subtree a is entirely within subtree b and make sure that a and b are not the same
        ))
      )
    ))

  subtree_distances.filtered.not_nested <- subtree_distances.filtered %>%
    filter(nested == FALSE)

  if (plot_clades_to_prune == TRUE) {
    #plot this on the tree to test
    temp_tree <- ggtree::ggtree(input_tree, right = T, ladderize = T)

    for (node.num in subtree_distances.filtered.not_nested$subtree_number) {
      temp_tree <- temp_tree +
        ggtree::geom_cladelabel(node = as.numeric(node.num) + length(input_tree$tip.label),
                                label = as.factor(node.num))
    }

    #write out the tree
    ggplot2::ggsave(plot = temp_tree,
                    filename = "original_tree_with_clades_to_prune.pdf",
                    height = 12,
                    width = 8)
  }



  #6 for all remaining, select one leaf per subtree, and remove the others
  selected_tips <- subtree_distances.filtered.not_nested %>%
    #select one tip per subtree
    mutate(selected = purrr::map(tips, ~ sample(.x, 1))) %>%
    # create list of all the rest
    mutate(to_remove = purrr::map2(tips, selected, function(a,b) a[!a %in% b])) %>%
    pull(to_remove) %>%
    unlist()


  pruned_tree <- micro.gen.extra::remove_tips_from_tree(tree = input_tree,
                                                        tips = selected_tips)

  #7 return the altered tree
  return(pruned_tree)

}






