
#' Create a dnaseqs list of dataframes for genoplotR
#'
#' Create a dnaseqs list of dataframes for genoplotR.
#' Uses a gggenes dataframe outputted from micro.gen.extra::gggenes_df_from_gff_dir().
#' Colours for genes can be set.
#'
#' @param gggenes_df dataframe created from micro.gen.extra::gggenes_df_from_gff_dir()
#' @param set_gene_fill_colours TRUE/FALSE as to whether genes should be coloured. If TRUE, the variable name for the
#' colour must be stored in a column called 'colour_variable' in the gggenes_df [Default = FALSE]
#' @param gene_colours Provide a named vector of hexadecimal colours,
#' where names correspond to the values in gggenes_df$colour_variable [Default = FALSE]
#' @param gene_rarity Minimum gene frequency across entire dataframe required to colour the gene, provided as an integer
#' (all genes lower than this frequency will be coloured white) [Default = 1]
#' @param clean_up_files TRUE/FALSE as to whether intermediate working files are removed [Default = TRUE]
#'
#' @return A dnaseqs list of dataframes for genoplotR
#' @examples
#' #Basic:
#' create_dnaseqs(input_gggenes_df)
#'
#' #More advanced:
#' create_dnaseqs(input_gggenes_df, TRUE, provided_colours_vector, 5, FALSE)
#'
#' #This will create a dnaseqs list of dataframes where any genes found at least five times across the dataset
#' #will be coloured according the colouring scheme set in "provided_colours_vector". Files created during the
#' #process of the function will be kept.
#'



create_dnaseqs <- function(gggenes_df,set_gene_fill_colours = FALSE, gene_colours = "random", gene_rarity = 1,clean_up_files = TRUE) {


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



  ##set some fills if requested:

  if (set_gene_fill_colours == T) {

    #check for colour_variable in colnames:


    if (!("colour_variable" %in% colnames(gggenes_df))) { # check if there is a column for the colour variable
      stop("Please have a column named 'colour_variable' in the 'gggenes_df' input")
    }

    #now create random colours if this is selected

    if (gene_colours[1] == "random") {

      gene_colours_df <- data.frame(table(as.character(gggenes_df$colour_variable))) %>%
        rename(colour_variable = Var1, occurence = Freq) %>%
        mutate(colour = random_n_colours(length(colour_variable),T))%>%
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

  for (gff_file in unique(gggenes_df$filename_prefix)) {
    temp_df <- subset.data.frame(gggenes_df, filename_prefix == gff_file)


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
    gff_file <- as.character(unique(value$filename_prefix))
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
    dna_seq_objects[[gff_file]] <- genoPlotR::read_dna_seg_from_tab(paste0(gff_file,".tab"), header = T)
    ### fix colours - had to remove hash for it to easily readable:

    #check that the colour is hex, if so - add a hash back onto the front of the hex number, otherwise leave as the string, e.g. "red", "blue", "pink" etc
    dna_seq_objects[[gff_file]]$fill <- purrr::map(dna_seq_objects[[gff_file]]$fill, ~ ifelse(is.hex(.x, hash = FALSE), paste('#', .x , sep=""), .x))

    #fix NA:
    dna_seq_objects[[gff_file]]$fill <- gsub("#NA",NA,dna_seq_objects[[gff_file]]$fill)
    ### and fix the colours for the outline of the genes:
    #if the colour is set as black as a hex (i.e. 000000), then read_dna_seg_from_tab will read it in as just 0 - therefore there is no colour
    dna_seq_objects[[gff_file]]$col <- purrr::map(dna_seq_objects[[gff_file]]$col, ~ ifelse(.x == 0, "000000", .x)) # therefore, if the colour is just 0, then change this to "000000"
    #check that the colour is hex, if so - add a hash back onto the front of the hex number, otherwise leave as the string, e.g. "red", "blue", "pink" etc
    dna_seq_objects[[gff_file]]$col <- purrr::map(dna_seq_objects[[gff_file]]$col, ~ ifelse(is.hex(.x, hash = FALSE), paste('#', .x , sep=""), .x))

    dna_seq_objects[[gff_file]]$col <- gsub("#NA",NA,dna_seq_objects[[gff_file]]$col)

    #clean up if requested:
    if (clean_up_files == T) {
      file.remove(paste0(gff_file,".tab"))
    }

  }


  return(dna_seq_objects)

}



#' Create blast comparisons for synteny plots
#'
#' Create blast comparisons for synteny plots provided a list of fasta sequences,
#' These can then be used as a readable format for genoplotR for synteny plots.
#' The order of the fasta sequences will be the order of the genomic regions to be compared in the synteny plot
#' from top to bottom.
#' blastn and makeblastdb
#'
#' @param blast_order List of the fasta files to be compared, not including the ".fasta" extension.
#' @param fasta_dir Directory where the fasta files are kept. If in the working directory, set this to "."
#' @param blast_results_dir Directory to put the blast results into [Default = "blast_results"]
#' @param overwrite TRUE/FALSE to overwrite a previously existing blast_results_dir [Default = FALSE]
#' @param color_scheme_name color scheme to be used for the genoplotR comparison [Default = "red_blue"]
#'
#' @return a dna_comparisons list for genoplotR in addition to a directory with blast comparisons for a list of fasta sequences.
#' Comparisons will be sequential. e.g. sequence_1 vs sequence_2, sequence_2 vs sequence_3, sequence_3 vs sequence_4
#'
#' @examples
#' create_blast_comparisons(blast_order_list, "dir_with_fasta_files")


create_blast_comparisons <- function(blast_order,fasta_dir,blast_results_dir = "blast_results", overwrite = FALSE, color_scheme_name = "red_blue") {

  comparisons <- list() # empty list for comparisons

  #make directory for blast_files
  if (dir.exists(blast_results_dir)) {
    if (isTRUE(overwrite)) {
      unlink(blast_results_dir)
    } else {
      stop(glue::glue("{blast_results_dir} is already a directory. Exiting."))
    }
  } else {
    dir.create(blast_results_dir)
  }

  for (i in 1:(length(blast_order)-1)) { # loop through and make pairwise blast comparisons between n and n+1
    #cluster <- names(clusters)[i]

    reference <- blast_order[i] #gsub(".gff","",clusters[[i]])
    query <- blast_order[i+1] #gsub(".gff","",clusters[[i+1]])

    system(glue::glue("makeblastdb -in {fasta_dir}/{query}.fasta -dbtype 'nucl' -title tmp_database -out tmp_database -parse_seqids"))
    system(glue::glue("blastn -task blastn -db tmp_database -perc_identity 20 -query {fasta_dir}/{reference}.fasta -evalue 10000 -outfmt 6 -out {blast_results_dir}/{query}_vs_{reference}"))

    #### see if getting rid of hashes works:
    system(glue::glue("sed -i 's/#/_/g' {blast_results_dir}/{query}_vs_{reference}"))

    ### load comparison into genoplotR data
    comparisons[[i]] <- try(genoPlotR::read_comparison_from_blast(glue::glue("{blast_results_dir}/{query}_vs_{reference}")))
    comparisons[[i]] <- subset.data.frame(comparisons[[i]], aln_len > 100 & e_value < 10)
  }


  #6 plot?
  for ( i in seq(1,length(comparisons))) {
    comparisons[[i]]$col <- genoPlotR::apply_color_scheme(comparisons[[i]]$per_id,
                                               direction=comparisons[[i]]$direction,
                                               color_scheme=color_scheme_name,
                                               rng=c(30,100))
  }

  #clean up the tmp_database_files:

  for (filename in list.files(path = ".", pattern = "tmp_database")) {
    file.remove(filename)
  }


  return(comparisons)

}


#' Plot the key for a genoplotR synteny plot
#'
#' Plot the key for a genoplotR synteny plot, output as svg for both forward and reverse blast hits.
#' Currently only a red/blue color scheme available
#'
#' @param color_scheme_name Colour scheme to be used. [Default = "red_blue"]
#'
#' @return Svg plots with individual keys for forward and reverse blast hits
#' @examples
#' plot_genoplotR_key()
#'
#' plot_genoplotR_key(color_scheme_name = "red_blue")

plot_genoplotR_key <- function(color_scheme_name = "red_blue") {




  #first need to add the color.bar function
  color.bar <- function(lut, min, max, axis_bool, alpha_num, ylab_string, title_string, nticks=8, ticks=seq(min, max, len=nticks)) {
    scale = (length(lut)-1)/(max-min)
    #dev.new(width=1.75, height=5)
    if (axis_bool == T) {
      plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', yaxt='n', xlab='', main = title_string, ylab=ylab_string)
      axis(2, ticks, las=1)
    } else {
      plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', yaxt='n', xlab='', ylab='', main = title_string)
    }
    for (i in 1:(length(lut)-1)) {
      y = (i-1)/scale + min
      rect(0,y,10,y+1/scale, col=alpha(lut[i], alpha_num), border=NA)
    }
  }




  if (color_scheme_name %in% c("red_blue")) {
    invisible()
  } else {
    stop(glue::glue("{color_scheme_name} not in available color schemes Exiting."))
  }




  ##### make scales:
  colour_nums <- seq(30,100, by=1)
  colours_act_forward <- genoPlotR::apply_color_scheme(colour_nums,
                                            direction=1,
                                            color_scheme="red_blue",
                                            rng=c(30,100))
  colours_act_reverse <- genoPlotR::apply_color_scheme(colour_nums,
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


#' Plot a genoplotR synteny plot
#'
#' Plot a genoplotR synteny plot, using the dnaseqs from micro_gen_extra::create_dnaseqs()
#' and the blast comparisons from micro_gen_extra::create_blast_comparisons()
#'
#'
#' @param input_genoplotR_set Input list of the dna_seqs and the blast comparisons
#' @param output_prefix Prefix for the output files (svg and png)
#'
#' @return A png and svg file of synteny plots produced using genoplotR
#' @examples
#' plot_genoplot_data(input_genoplotR_data, "output_prefix_name")
#'
#' More specifically:
#'
#' #First create an empty list
#' genoplotr_data <- list()
#'
#' #Add the dnaseqs:
#' genoplotr_data$dnaseqs <- micro.gen.extra::create_dnaseqs(input_gggenes_df)
#'
#' #Then add the blast comparisons:
#' genoplotr_data$comparisons <- micro.gen.extra::create_blast_comparisons(blast_order_list, "dir_with_fasta_files")
#'
#' #Then plot:
#' plot_genoplot_data(genoplotr_data, "output_prefix_name")
#'
#' #To add colours, and change blast parameters (to be added), see the documentation for micro_gen_extra::create_dnaseqs
#' #and micro_gen_extra::create_blast_comparisons




plot_genoplot_data <- function(input_genoplotR_set,output_prefix) {

  max_length <- 0
  for (c in input_genoplotR_set$dna_seqs) {
    if (max(c$end) > max_length) {
      max_length <- max(c$end)
    }
  }


  #set the heights

  num_comparisons <- length(input_genoplotR_set$dna_seqs)

  output_width <- (max_length / 10000 ) * 3

  output_height <- (num_comparisons ^ 2) / 4

  #output <- "test_output"

  png(glue::glue("{output_prefix}.png"), res = 300, units = "in", width = output_width, height = output_height,pointsize = 8)
  genoPlotR::plot_gene_map(input_genoplotR_set$dna_seqs, input_genoplotR_set$comparisons,
                scale=T,
                legend = T,
                arrow_head_len = 400,
                annotation_height = 0.5)
  dev.off()

  svg(glue::glue("{output_prefix}.svg"), width = output_width, height = output_height,pointsize = 8)
  genoPlotR::plot_gene_map(input_genoplotR_set$dna_seqs, input_genoplotR_set$comparisons,
                scale=T,
                legend = T,
                arrow_head_len = 400,
                annotation_height = 0.5)
  dev.off()
}


#' Take a list of "x" gff3 files and create a synteny plot
#'
#' Take gff files, read annotation and extract fasta sequence.
#' Then create genoplotR data structures, then plot genoplotR synteny plot.
#' This uses the dnaseqs from micro_gen_extra::create_dnaseqs()
#' and the blast comparisons from micro_gen_extra::create_blast_comparisons()
#'
#'
#' @param gff_list Input list of gff files (with extensions to paths)
#' @param output_prefix Prefix for the output files (svg and png)
#' @param temp_dir Directory to put fasta files into [Default = "temp_fasta"]
#' @param clean_up TRUE/FALSE - clean up the temp fasta directory [Default = TRUE]
#'
#' @return A png and svg file of synteny plots produced using genoplotR
#' @examples
#' Insert examples here




plot_synteny_gff_files <- function(gff_list,output_prefix = "synteny_plot", temp_dir = "temp_fasta", clean_up = TRUE) {


  names(gff_list) <- unlist(purrr::map(gff_list, ~ tail(unlist(stringr::str_split(.x,"/")),1)))


  gggenes_df <- micro.gen.extra::gggenes_df_from_gff_list(gff_list = gff_list)


  genoplotR_data <- list()

  #now get the dna_seqs from the gggenes dataframe
  genoplotR_data$dna_seqs <- micro.gen.extra::create_dnaseqs(gggenes_df,set_gene_fill_colours = FALSE)

  micro.gen.extra::fasta_from_gff_list(gff_list = gff_list,
                                       gff_names = names(gff_list),
                                       fasta_dir = temp_dir,
                                       clean_up = FALSE)

  #run the blast comparisons and read them for genoplotR
  genoplotR_data$comparisons <- micro.gen.extra::create_blast_comparisons(names(genoplotR_data$dna_seqs), temp_dir,overwrite = TRUE)

  svg("synteny_selected.svg",width = 15, height = 30)
  genoPlotR::plot_gene_map(genoplotR_data$dna_seqs, genoplotR_data$comparisons,
                           scale=T,
                           legend = T,
                           arrow_head_len = 400,
                           annotation_height = 0.5,
                           plot_new = FALSE)
  dev.off()

  #clean up
  if (isTRUE(clean_up)) {
    unlink(temp_dir, recursive = TRUE)
  }

  return(genoplotR_data)

}
