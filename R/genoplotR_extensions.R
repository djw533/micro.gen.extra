create_dnaseqs <- function(gggenes_df,set_gene_fill_colours,gene_colours,gene_rarity,clean_up_files) {


  # where :
  # gggenes_df == a dataframe of all genes across all gene clusters to be compared as created by gggenes_df_from_gff_dir
  # set_gene_fill_colours == TRUE or FALSE for whether you want a fill to be set for the genes. The variable for the colour should
  #   be in the gggenes_df, in a column named colour_variable (n.b this is the variable to be used for colouring, specifically not the colour
  #   to be used! If set to FALSE, all genes will be white.
  # gene_colours == Either a named vector or FALSE.  Set this to either a named vector with hex colours.
  #   e.g. c("variable_1" = "#BEBEBE", "variable_2" = "#AE473E"). If you want random colours for each of the variables - set this to FALSE
  # gene_rarity == A number which will not colour genes that appear less than that value across all gene clusters. e.g. if this is set to 5,
  #   there are 10 different gene clusters, and a certain variable (set in the column colour_variable) appears only 4 times across all 10
  #   gene clusters, then this gene will just be coloured white regardless.  This is set so that if you want to focus on genes that are
  #   shared across all gene clusters, then very rare genes won't be highlighted. (The effect from this will probably also be accentuated if
  #   creating a synteny plot, as this gene will likely not be shared with gene clusters above and below).

  # gggenes_df = gggenes_w_cdhit_clusters
  # clustered_genes = TRUE
  # gene_rarity = 5
  #


  ##set some fills if requested:

  if (set_gene_fill_colours == T) {


    #now create random colours if this is selected

    if (isFALSE(gene_colours) == F) {

      gene_colours_df <- data.frame(table(as.character(gggenes_df$colour_variable))) %>%
        rename(colour_variable = Var1, occurence = Freq) %>%
        mutate(colour = random_x_colours(length(colour_variable),T))%>%
        #now set all colours top grey if less than gene rarity:
        mutate(colour = ifelse(occurence < gene_rarity, "#BEBEBE", colour))


      ### now put these into a col vector:
      gene_colours <- gene_colours_df$colour
      names(gene_colours) <- gene_colours_df$colour_variable
    }

    # no else required - we set the gene_colours in the function call
    # but could perhaps check the contents of the gene colours list that has been produced here.
    # should be a named character vector, that has gene colours for all of the genes in the dataframe, and all should be hexadecimal numbers.

  }

  #set colours:
  feature_colours <- c("CDS" = "#000000" ,
                       "ncRNA" =  "#0000ff",
                       "tRNA" = "#ff0000")


  #split gggenes_df into a list of dataframes:

  all_gggenes_data <- list()

  for (gff_file in unique(gggenes_df$gff_filename)) {
    temp_df <- subset.data.frame(gggenes_df, gff_filename == gff_file)


    all_gggenes_data[[gff_file]] <- temp_df
  }


  dna_seqs <- list()
  dna_seq_objects <- list()

  ###now write into the dna_seqs objects:
  for (i in 1:length(all_gggenes_data)) {
    #set names again
    #gff_file <- gsub(".gff","",as.character(clusters[i])) #  changing this below so that it is the cluster number!:
    value <- all_gggenes_data[[i]]
    ##get the filename of the gff in question:
    gff_file <- as.character(unique(value$gff_filename))
    ##put data into the dnaseqs list of dfs
    dna_seqs[[gff_file]] <- data.frame(matrix(ncol = 18, nrow = nrow(value)))
    colnames(dna_seqs[[gff_file]]) <- c("name","start","end","strand","length","pid","gene","synonym","product","proteinid","feature","gene_type","col","fill","lty","lwd","pch","cex")
    ### put variables into this df from the value :
    dna_seqs[[gff_file]]$name <- gsub("#","_",value$gene)
    dna_seqs[[gff_file]]$start <- value$start
    dna_seqs[[gff_file]]$end <- value$end
    dna_seqs[[gff_file]]$strand <- value$direction
    dna_seqs[[gff_file]]$length <- (value$end - (value$start - 1)) / 3
    dna_seqs[[gff_file]]$pid <- NA
    dna_seqs[[gff_file]]$gene <- gsub("#","_",value$gene) # gsub out any hash
    dna_seqs[[gff_file]]$synonym <- NA
    dna_seqs[[gff_file]]$product <- NA
    dna_seqs[[gff_file]]$proteinid <- NA
    dna_seqs[[gff_file]]$feature <- value$type
    dna_seqs[[gff_file]]$gene_type <- "arrows"
    dna_seqs[[gff_file]]$lty <- 1
    dna_seqs[[gff_file]]$lwd <- 1
    dna_seqs[[gff_file]]$pch <- 8
    dna_seqs[[gff_file]]$cex <- 1

    ##set colours if requested:
    if (set_gene_fill_colours == T) {
      dna_seqs[[gff_file]]$col <- gsub("#","", unname(feature_colours[match(value$type, names(feature_colours))]))      # change colour depending on the type of the feature??
      dna_seqs[[gff_file]]$fill <- gsub("#","", unname(gene_colours[match(value$colour_variable, names(gene_colours))]))
    } else {
      dna_seqs[[gff_file]]$col <- gsub("#","", unname(feature_colours[match(value$type, names(feature_colours))]))      # change colour depending on the type of the feature??
      dna_seqs[[gff_file]]$fill <- "FFFFFF"
    }


    write.table(file = paste0(gff_file,".tab"), dna_seqs[[gff_file]] , sep = "\t", quote = F, row.names = F)
    dna_seq_objects[[gff_file]] <- read_dna_seg_from_tab(paste0(gff_file,".tab"), header = T)
    ### fix colours - had to remove hash for it to easily readable:
    dna_seq_objects[[gff_file]]$fill <- paste('#',dna_seq_objects[[gff_file]]$fill, sep="")
    #fix NA:
    dna_seq_objects[[gff_file]]$fill <- gsub("#NA",NA,dna_seq_objects[[gff_file]]$fill)
    ### and fix the colours for the outline of the genes:
    #if the colour is set as black as a hex (i.e. 000000), then read_dna_seg_from_tab will read it in as just 0 - therefore there is no colour
    dna_seq_objects[[gff_file]]$col <- map(dna_seq_objects[[gff_file]]$col, ~ ifelse(.x == 0, "000000", .x)) # therefore, if the colour is just 0, then change this to "000000"
    dna_seq_objects[[gff_file]]$col <- paste('#',dna_seq_objects[[gff_file]]$col, sep="") # then paste the hash back in
    dna_seq_objects[[gff_file]]$col <- gsub("#NA",NA,dna_seq_objects[[gff_file]]$col)

    #clean up if requested:
    if (clean_up_files == T) {
      file.remove(paste0(gff_file,".tab"))
    }

  }


  return(dna_seq_objects)

}


create_blast_comparisons <- function(blast_order,fasta_dir,blast_results_dir) {

  comparisons <- list() # empty list for comparisons

  #make directory for blast_files
  if (dir.exists(blast_results_dir)) {
    stop(glue("{blast_results_dir} is already a directory. Exiting."))
  } else {
    dir.create(blast_results_dir)
  }

  for (i in 1:(length(blast_order)-1)) { # loop through and make pairwise blast comparisons between n and n+1
    #cluster <- names(clusters)[i]

    reference <- blast_order[i] #gsub(".gff","",clusters[[i]])
    query <- blast_order[i+1] #gsub(".gff","",clusters[[i+1]])

    system(glue("makeblastdb -in {fasta_dir}/{query}.fasta -dbtype 'nucl' -title tmp_database -out tmp_database -parse_seqids"))
    system(glue("blastn -task blastn -db tmp_database -perc_identity 20 -query {fasta_dir}/{reference}.fasta -evalue 10000 -outfmt 6 -out {blast_results_dir}/{query}_vs_{reference}"))

    #### see if getting rid of hashes works:
    system(glue("sed -i 's/#/_/g' {blast_results_dir}/{query}_vs_{reference}"))

    ### load comparison into genoplotR data
    comparisons[[i]] <- try(read_comparison_from_blast(glue("{blast_results_dir}/{query}_vs_{reference}")))
    comparisons[[i]] <- subset.data.frame(comparisons[[i]], aln_len > 100 & e_value < 10)
  }


  #6 plot?
  for ( i in seq(1,length(comparisons))) {
    comparisons[[i]]$col <- apply_color_scheme(comparisons[[i]]$per_id,
                                               direction=comparisons[[i]]$direction,
                                               color_scheme="red_blue",
                                               rng=c(30,100))
  }

  #clean up the tmp_database_files:

  for (filename in list.files(path = ".", pattern = "tmp_database")) {
    file.remove(filename)
  }


  return(comparisons)

}



plot_genoplotR_key <- function(color_scheme) {

  if (color_scheme %in% c("red_blue")) {
    invisible()
  } else {
    stop(glue("{color_scheme} not in available color schemes Exiting."))
  }




  ##### make scales:
  colour_nums <- seq(30,100, by=1)
  colours_act_forward <- apply_color_scheme(colour_nums,
                                            direction=1,
                                            color_scheme="red_blue",
                                            rng=c(30,100))
  colours_act_reverse <- apply_color_scheme(colour_nums,
                                            direction=-1,
                                            color_scheme="red_blue",
                                            rng=c(30,100))


  ### remove the alpha:
  colours_act_forward <- gsub("80","",colours_act_forward)
  colours_act_reverse <- gsub("80","",colours_act_reverse)

  svg("forward_scale.svg")#, res = 300, units = "cm", width = 4, height = 8,pointsize = 8)
  forward_ramp <- color.bar((colorRampPalette(colours_act_forward))(length(colour_nums)), 30, 100, T, 0.5, "Percentage identity (%)", "Forward")
  dev.off()

  svg("reverse_scale.svg")#, res = 300, units = "cm", width = 4, height = 8,pointsize = 8)
  reverse_ramp <- color.bar((colorRampPalette(colours_act_reverse))(length(colour_nums)), 30, 100, F, 0.5, "Percentage identity (%)", "Reverse")
  dev.off()



}


plot_genoplot_data <- function(input_genoplotR_set,output_prefix) {

  max_length <- 0
  for (c in genoplotr_data$dna_seqs) {
    if (max(c$end) > max_length) {
      max_length <- max(c$end)
    }
  }


  #set the heights

  num_comparisons <- length(genoplotr_data$dna_seqs)

  output_width <- (max_length / 10000 ) * 3

  output_height <- (num_comparisons ^ 2) / 4

  #output <- "test_output"

  png(glue("{output_prefix}.png"), res = 300, units = "in", width = output_width, height = output_height,pointsize = 8)
  plot_gene_map(input_genoplotR_set$dna_seqs, input_genoplotR_set$comparisons,
                scale=T,
                legend = T,
                arrow_head_len = 400,
                annotation_height = 0.5)
  dev.off()

  svg(glue("{output_prefix}.svg"), width = output_width, height = output_height,pointsize = 8)
  plot_gene_map(input_genoplotR_set$dna_seqs, input_genoplotR_set$comparisons,
                scale=T,
                legend = T,
                arrow_head_len = 400,
                annotation_height = 0.5)
  dev.off()
}


