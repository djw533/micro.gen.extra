get_intersect_members <- function (upset_df, groups, exclusive = TRUE){


  # try to work out the rowsums

  all_rowssums <- upset_df %>%
    rowwise %>%
    mutate(row.sum = sum(c_across(where(is.numeric))))


  if (isTRUE(exclusive)) {

    selected_rowssums <- upset_df %>%
      tidyr::pivot_longer(!Gene, names_to = "group", values_to = "presence") %>%
      filter(group %in% groups) %>%
      tidyr::pivot_wider(names_from = "group", values_from = "presence") %>%
      rowwise %>%
      mutate(row.sum.filtered = sum(c_across(where(is.numeric)))) %>%
      full_join(all_rowssums) %>%
      mutate(intersection = ifelse( (row.sum.filtered == row.sum & row.sum.filtered > 0 & row.sum.filtered == length(groups)), T, F)) %>% # find instances where the rowsums for the selected groups and the rowsums for all groups are the same and also more than 0
      filter(isTRUE(intersection))
  }

  if ( isFALSE(exclusive)) {
    selected_rowssums <- upset_df %>%
      tidyr::pivot_longer(!Gene, names_to = "group", values_to = "presence") %>%
      filter(group %in% groups) %>%
      tidyr::pivot_wider(names_from = "group", values_from = "presence") %>%
      rowwise %>%
      mutate(row.sum.filtered = sum(c_across(where(is.numeric)))) %>%
      full_join(all_rowssums) %>%
      mutate(intersection = ifelse( row.sum.filtered == length(groups), T, F)) %>% # find instances where the rowsums for the selected groups and the rowsums for all groups are the same and also more than 0
      filter(isTRUE(intersection))
  }

  return(selected_rowssums$Gene)

}





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


