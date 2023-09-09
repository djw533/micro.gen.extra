#' Get members (genes) of a set intersection shown in an UpsetR plot
#'
#' Given a list of groups in a set analysis plotted in UpsetR, get all members
#' (i.e. Genes in this instance) in the intersection of those sets. Built to 
#' extract genes in a panaroo dataset.
#'
#' @param upset_df  UpsetR dataframe (with the members column labelled as 
#'                  "Gene")
#' @param groups    List of the groups to get the intersection for
#' @param exclusive TRUE/FALSE. If TRUE, only get the intersection of ONLY the 
#'                  groups in groups.
#'                  If FALSE, return all values shared by the groups in groups,
#'                  irrespective of whether those values are also in other 
#'                  groups. [Default = TRUE]
#'
#' @return List of all members in the intersection
#' @examples
#' get_intersect_members(dataframe_used_for_upsetR_plot, list_of_groups)

get_intersect_members <- function (
    upset_df, groups, ids = "Gene", exclusive = TRUE
){

  # try to work out the rowsums
  all_rowssums <- upset_df %>%
    #rowwise %>%
    rowwise() %>%
    mutate(row.sum = sum(c_across(where(is.numeric))))
  
  
  upset_df %>%
    tidyr::pivot_longer(
      cols      = !all_of(ids),
      names_to  = "group",
      values_to = "presence"
    ) %>%
    filter(group %in% groups) %>%
    tidyr::pivot_wider(names_from = "group", values_from = "presence") %>%
    rowwise() %>%
    mutate(row.sum.filtered = sum(c_across(where(is.numeric)))) %>%
    full_join(all_rowssums) %>%
    # Find instances where the rowsums for the selected groups and the rowsums 
    # for all groups are the same and also more than 0
    mutate(intersection = case_when(
      isTRUE(exclusive)  ~ ifelse(
        test = row.sum.filtered == row.sum &
          row.sum.filtered > 0 &
          row.sum.filtered == length(groups),
        yes  = T,
        no   = F
      ),
      isFALSE(exclusive) ~ ifelse(
        test = row.sum.filtered == length(groups),
        yes  = T,
        no   = F
      )
    )) %>%
    filter(isTRUE(intersection)) %>% 
    pull(ids)
}



#' Plot gene accumulation curves
#'
#' Plot gene accumulation curves from a gene_presence_absence file from a pan-genome analysis
#' using vegan and ggplot2
#'
#' @param gene_presence_absence_file gene_presence_absence.Rtab file from Roary or Panaroo
#' @param groups Named vector with a list of lists of strains, with the name being each cluster / lineage / species etc.
#'
#' @return Plot of gene accumulation curves, and the dataframe used to plot these
#' @examples
#' plot_gene_accumulation_curves(gene_presence_absence.Rtab, list_of_groups)


plot_gene_accumulation_curves <- function(gene_presence_absence_file,groups) {



  pangenome_df <- read.table(gene_presence_absence_file,
                             header = T,
                             sep="\t",
                             quote = "",
                             comment.char = "",
                             stringsAsFactors = F,
                             check.names = F)

  #get the strains in the pangenome:
  pangenome_strains <- data.frame("strain" = colnames(pangenome_df)[2:(ncol(pangenome_df))])

  # join the groups into the pangenome:

  pangenome_strains_w_group <- pangenome_strains %>%
    left_join(groups, by = "strain") %>%
    filter(! is.na(group))



  ### print out warning if some of the pangenome strains have now been removed from the analysis as they aren't included in the groups df:

  if (nrow(pangenome_strains_w_group) < nrow(pangenome_strains)) {
    print(glue::glue("Warning - {as.character(nrow(pangenome_strains) - nrow(pangenome_strains_w_group))} removed from analysis as strains were not present in the groups"))

    # get the removed strains as a list:
    removed_strains <- pangenome_strains %>%
      left_join(groups, by = "strain") %>%
      filter(is.na(group)) %>%
      select("strain")

    print(glue::glue("Removed {removed_strains$strain}"))
  }

  #use vegan to create the pangenome accumulation curves:
  df_phylogroups = data.frame(cluster = character(0),
                              richness = numeric(0),
                              genomes = numeric(0),
                              sd = numeric(0))


  for (group_name in unique(pangenome_strains_w_group$group)){
    print(group_name)
    ## get strains for this phylogroup:
    group_strains <- pangenome_strains_w_group %>%
      filter(group == group_name)

    #subset out the diff gene_presence_absence_data:
    temp_df <- pangenome_df %>%
      select(group_strains$strain)

    ## set the rownames as genes:
    rownames(temp_df) <- pangenome_df$Gene


    if (ncol(temp_df) > 2 ) { # this will only work if more than 2 strains in the poppunk cluster
      no_values <- data.frame(t(temp_df[apply(temp_df[,-1], 1, function(x) !all(x==0)),]))
    } else {
      next
    }

    ##add to the accum dataset:
    sp <- vegan::specaccum(no_values, "random", permutations=100)

    df_phylogroups = rbind(df_phylogroups, data.frame(cluster = rep(group_name, length(sp$sites)),
                                                      richness = sp$richness,
                                                      genomes = sp$sites,
                                                      sd = sp$sd))

  }


  #clean up and add min and max:
  df_phylogroups = cbind(df_phylogroups, min= df_phylogroups$richness-df_phylogroups$sd, max = df_phylogroups$richness+df_phylogroups$sd )
  df_phylogroups = df_phylogroups[which(df_phylogroups$genomes %% 1 == 0),] ## to visualise more clearly


  ggplot(df_phylogroups, aes(x = genomes, y = richness, color = as.factor(cluster)))+ geom_line(alpha = 0.6, size = 1) +#+  geom_point(size = 2, alpha = 0.8) +
    theme_bw(base_size = 16) + xlab("Number of genomes") +
    ylab("Unique CDSs") +
    geom_ribbon(aes(ymin = min, ymax = max, fill = as.factor(cluster)), alpha  = 0.1, color = NA) +
    #facet_wrap(vars(phylogroup)) +
    #scale_color_manual(name = "Species",values = c(tree_cols)) +
    #  scale_fill_manual(name = "Species",values = c(tree_cols)) +#+ scale_shape_manual(values = graphics$Shape, guide = F)  +
    theme(legend.position = "right")


  return(df_phylogroups)

}





#' Create UpSetR dataframes from a gene_presence_absence.Rtab file
#'
#' Create UpSetR dataframes from a gene_presence_absence.Rtab file, with the different groups to create the
#' sets being defined by clusters - e.g. phylogroups/lineages calculated by fastANI/PopPUNK/fastbaps etc
#'
#' @param gene_presence_absence_file gene_presence_absence.Rtab file from Roary or Panaroo
#' @param groups Dataframe of the groups with two columns, one named "strain" listing strain names, and the other named "group" listing the group/lineage etc that each strain belongs to
#' @param levels Levels to define the pan-genome [Default = c(0,0.15,0.95,1)]
#' @param labels Labels for groups created between each of the values in levels. Should be of length : (length(levels) - 1).
#' [Default = c("Cloud","Shell","Core")]
#' @return UpsetR dataframe ready for plotting
#' @examples
#' create_upset_dfs(gene_presence_absence.Rtab, list_of_groups)

create_upset_dfs <- function(gene_presence_absence_file,groups,levels = c(0,0.15,0.95,1), labels = c("Cloud","Shell","Core")) {

 #function to get position of the gene within the given levels

  position_number <- function(x, levels, labels) {

    #and for x:
    if ( ! is.numeric(x)) {
      stop("x must be numeric")
    }

    #function

    for (i in 1:(length(levels)-1)) {

      min = levels[i]
      max = levels[i+1]

      if (x == 0) {
        label = "Absent"
        return(label)
      }

      if (x >= min & x < max) {
        label = labels[i]

        return(label)
      }

      if ( x == levels[length(levels)]) { # if x is at the top limits of the levels - then but in the top level ( fix for above x < max check)
        label = labels[length(labels)]
        return(label)
      }

    }

  }

  ###

  #check length of levels is 1 less than the length of the labels:

  if ( (length(levels) - length(labels)) != 1) {
    stop("Levels length should be 1 greater than labels")
  }

  #check all levels are between 0 and 1:

  for (level_num in levels) {
    if (level_num < 0 | level_num > 1) {
      stop("Levels must all be within 0-1 (inclusive")
    }
  }

  #finally check that the levels increase in number:

  for (i in 1:(length(levels)-1) ) {
    if (levels[i] >= levels[i+1]) {
      stop("Numbers in levels must ascend")
    }
  }

  pangenome_df <- read.table(gene_presence_absence_file,
                             header = T,
                             sep="\t",
                             quote = "",
                             comment.char = "",
                             stringsAsFactors = F,
                             check.names = F)


  ###try using a cleaner dplyr method?

  gene_occurence_df <- pangenome_df %>%
    tidyr::pivot_longer(!Gene,names_to="strain",values_to = "presence") %>%
    left_join(groups, by = "strain") %>% # add groups in
    group_by(group,Gene) %>%
    summarise(perc_in_group = sum(presence) / length(strain)) %>% # get incidence of the gene within each lineage
    #filter(perc_in_group > 0) %>% # remove genes not seen in individual groups of strains
    mutate(occurence_label = purrr::map(perc_in_group, ~ position_number(.x, levels, labels))) ## map over incidence and assign group for the gene


  #create the upset plot for each label:

  upset_dfs <- vector(mode = "list", length = length(labels))
  names(upset_dfs) <- labels

  for (i in 1:length(labels)) {

    l = labels[i]

    upset_df <- gene_occurence_df %>%
      mutate(perc_in_group = ifelse(occurence_label == l, 1, 0 )) %>%
      select(Gene,group,perc_in_group) %>%
      tidyr::pivot_wider(names_from = group, values_from = perc_in_group)

    #change the colnames back:
    colnames(upset_df) <- gsub(glue::glue("{l}[.]"),"",colnames(upset_df))

    upset_dfs[i] <-list(as.data.frame(upset_df))

  }

  return(upset_dfs)


}


plot_intersection_on_tree <- function(ggtree_object, intersection_members, upset_df, eggnog_annotation, pangenome_df_transposed, exclusive_int = TRUE, change_names = FALSE) {


  if (isTRUE(change_names)) {
    ##create annotation names
    annotation_names <- eggnog_annotation %>%
      filter(! "#" %in% Gene) %>%
      full_join(upset_df %>% select(Gene)) %>%
      mutate(heatmap_name = ifelse(is.na(desc), Gene, desc)) %>%
      #now arrange by the heatmap name - and add suffixes if the name is duplicated
      arrange(heatmap_name) %>%
      mutate(to_del = 1) %>%
      group_by(heatmap_name) %>%
      mutate(csum=cumsum(to_del)) %>%
      mutate(csum = ifelse(csum  == 1, "", paste0("_",csum,sep=""))) %>%
      mutate(heatmap_name2 = paste0(heatmap_name,csum, sep = "")) %>%
      pull(heatmap_name2, name = Gene)

  }

  # extract stuff
  intersection <- micro.gen.extra::get_intersect_members(upset_df = upset_df,  intersection_members, exclusive = exclusive_int)

  intersection_gheatmap <- pangenome_df_transposed %>%
    select(all_of(intersection))

  if (isTRUE(change_names)) {
    colnames(intersection_gheatmap) <- unname(annotation_names[c(colnames(intersection_gheatmap))])
  }

  #plot:
  ggtree_object %>% ggtree::gheatmap(intersection_gheatmap, color = NULL,
                                     colnames_position = "top",
                                     colnames_angle = 45,
                                     colnames_offset_y = 5,
                                     hjust = 0,
                                     width = 3,
                                     offset = 0.1) +
    scale_fill_gradient(low = "white", high = "blue") +
    theme(legend.position = "none",
          text = element_text(size = 2)) +
    ylim(0,max(t2$data$y)+150)

  return(intersection_gheatmap)

}



