#' Plot T6SS hamburger results for each strain
#'
#' Plot T6SS hamburger results for each strain
#'
#' @param hamburger_dir Hamburger output directory
#' @param strains List of strains to plot gene clusters from
#' @param output_dir Output directory to put image files
#' @param overwrite TRUE/FALSE for whether to overwrite output_dir (if it already exists)
#' @param format png/svg output, can select both - c("png","svg","both")
#' @param height Plot height in inches - default = 8
#' @param width Plot width in inches - default = 8
#'
#' @return png/svg files of gene arrow plots for extracted T6SS gene clusters using hamburger
#' @examples
#' gggenes_df_from_gff_dir("path/to/gff/directory")



plot_t6_clusters <- function(hamburger_dir, strains = "all", output_dir = ".", overwrite = F, format = "png", height = 8, width = 8) {

  #make output directory:

  if (! output_dir == ".") {
    if (dir.exists(output_dir)) {
      if (isTRUE(overwrite)) {
        unlink(output_dir, recursive = TRUE)
      } else {
        stop("Directory exists")
      }

    } else {
      dir.create(output_dir)
    }
  }

  # check output

  if (! format %in% c("svg","png","both")) {
    stop("Please choose one of 'png', 'svg', or 'both' for output")
  }

  #need some checks on the gggenes file here
  input <- read.csv(glue::glue("{hamburger_dir}/gggenes_input.csv"),
                    stringsAsFactors = F,
                    comment.char = "",
                    header = T)

  stats <- read.csv(glue::glue("{hamburger_dir}/cluster_stats.csv"),
                    stringsAsFactors = F,
                    comment.char = "",
                    header = T) %>%
    rename(operon = gene_cluster)

  ##add the contigs in
  input <- input %>%
    left_join(stats %>% select(operon,contig)) %>%
    mutate(new_cluster_name = paste(operon,contig,sep = "_"))


  if (! strains == "all") {
    input <- input %>%
      filter(strain %in% strains)
  }



  #set the colours:

  t6ss_cols <- c("TssA" = "#3cb44b",
                 "TssB" = "#ffe119",
                 "TssC" = "#e6194b",
                 "TssD" = "#4363d8",
                 "TssE" = "#ff1493",
                 "TssF" = "#911eb4",
                 "TssG" = "#46f0f0",
                 "TssH" = "#f032e6",
                 "TssI" = "#bcf60c",
                 "TssJ" = "#fabebe",
                 "TssK" = "#008080",
                 "TssL" = "#e6beff",
                 "TssM" = "#9a6324",
                 "PAAR_motif" = "#808000",
                 "Non-model" = "#ffffff")


  #plot gene clusters in each strain selected


  for (s in unique(input$strain)) {


    temp_input <- input %>%
      filter(strain == s)


    p <- ggplot2::ggplot(temp_input, ggplot2::aes(xmin = start, xmax = end,y = new_cluster_name, fill = gene, forward = direction), color = "black") +
      gggenes::geom_gene_arrow(arrowhead_height = grid::unit(3, "mm"),
                               arrow_body_height = grid::unit(2,"mm"),
                               arrowhead_width = grid::unit(2, "mm")) +
      #gggenes::geom_gene_label(aes(label = gene)) +
      ggplot2::geom_text(aes(x = (start + end) / 2,  label = CDS_identifier), size = 2, nudge_y = 0.1, angle = 45, hjust = 0) +
      scale_fill_manual(values=c(t6ss_cols)) +
      ggplot2::xlab(glue::glue("Position (nt)")) +
      ggplot2::ylab(glue::glue("Cluster name and contig")) +
      gggenes::theme_genes() +
      ggplot2::theme(legend.position = "bottom") +
      ggplot2::xlim(0, max(temp_input$end)*1.1)


    if (format == "both") {
      ggsave(plot = p,
             file = glue::glue("{output_dir}/{s}_T6SS_genes.png"),
             width = 8,
             height = 8,
             dpi = 600)
      ggsave(plot = p,
             file = glue::glue("{output_dir}/{s}_T6SS_genes.svg"),
             width = 8,
             height = 8,
             dpi = 600)


    } else {
      ggsave(plot = p,
             file = glue::glue("{output_dir}/{s}_T6SS_genes.{format}"),
             width = 8,
             height = 8,
             dpi = 600)

    }




  }


}




#' Plot hamburger results for each strain
#'
#' Plot T6SS hamburger results for each strain
#'
#' @param hamburger_dir Hamburger output directory
#' @param strains List of strains to plot gene clusters from
#' @param output_dir Output directory to put image files
#' @param overwrite TRUE/FALSE for whether to overwrite output_dir (if it already exists)
#' @param format png/svg output, can select both - c("png","svg","both")
#' @param height Plot height in inches - default = 8
#' @param width Plot width in inches - default = 8
#' @param colours Either "random" or a named vector list of colours for genes in the hmm query
#'
#' @return png/svg files of gene arrow plots for extracted gene clusters using hamburger
#' @examples
#' gggenes_df_from_gff_dir("path/to/gff/directory")



plot_hamburger <- function(hamburger_dir, strains = "all", output_dir = ".", overwrite = F, format = "png", height = 8, width = 8, colours = "random") {

  #make output directory:

  if (! output_dir == ".") {
    if (dir.exists(output_dir)) {
      if (isTRUE(overwrite)) {
        unlink(output_dir, recursive = TRUE)
      } else {
        stop("Directory exists")
      }

    } else {
      dir.create(output_dir)
    }
  }

  # check output

  if (! format %in% c("svg","png","both")) {
    stop("Please choose one of 'png', 'svg', or 'both' for output")
  }

  #need some checks on the gggenes file here
  input <- read.csv(glue::glue("{hamburger_dir}/gggenes_input.csv"),
                    stringsAsFactors = F,
                    comment.char = "",
                    header = T)

  stats <- read.csv(glue::glue("{hamburger_dir}/cluster_stats.csv"),
                    stringsAsFactors = F,
                    comment.char = "",
                    header = T) %>%
    rename(operon = gene_cluster)

  ##add the contigs in
  input <- input %>%
    left_join(stats %>% select(operon,contig)) %>%
    mutate(new_cluster_name = paste(operon,contig,sep = "_"))


  if (! strains == "all") {
    input <- input %>%
      filter(strain %in% strains)
  }



  #set the colours:

  t6ss_cols <- c("TssA" = "#3cb44b",
                 "TssB" = "#ffe119",
                 "TssC" = "#e6194b",
                 "TssD" = "#4363d8",
                 "TssE" = "#ff1493",
                 "TssF" = "#911eb4",
                 "TssG" = "#46f0f0",
                 "TssH" = "#f032e6",
                 "TssI" = "#bcf60c",
                 "TssJ" = "#fabebe",
                 "TssK" = "#008080",
                 "TssL" = "#e6beff",
                 "TssM" = "#9a6324",
                 "PAAR_motif" = "#808000",
                 "Non-model" = "#ffffff")


  #plot gene clusters in each strain selected


  for (s in unique(input$strain)) {


    temp_input <- input %>%
      filter(strain == s)


    p <- ggplot2::ggplot(temp_input, ggplot2::aes(xmin = start, xmax = end,y = new_cluster_name, fill = gene, forward = direction), color = "black") +
      gggenes::geom_gene_arrow(arrowhead_height = grid::unit(3, "mm"),
                               arrow_body_height = grid::unit(2,"mm"),
                               arrowhead_width = grid::unit(2, "mm")) +
      #gggenes::geom_gene_label(aes(label = gene)) +
      ggplot2::geom_text(aes(x = (start + end) / 2,  label = CDS_identifier), size = 2, nudge_y = 0.1, angle = 45, hjust = 0) +
      scale_fill_manual(values=c(t6ss_cols)) +
      ggplot2::xlab(glue::glue("Position (nt)")) +
      ggplot2::ylab(glue::glue("Cluster name and contig")) +
      gggenes::theme_genes() +
      ggplot2::theme(legend.position = "bottom") +
      ggplot2::xlim(0, max(temp_input$end)*1.1)


    if (format == "both") {
      ggsave(plot = p,
             file = glue::glue("{output_dir}/{s}_T6SS_genes.png"),
             width = 8,
             height = 8,
             dpi = 600)
      ggsave(plot = p,
             file = glue::glue("{output_dir}/{s}_T6SS_genes.svg"),
             width = 8,
             height = 8,
             dpi = 600)


    } else {
      ggsave(plot = p,
             file = glue::glue("{output_dir}/{s}_T6SS_genes.{format}"),
             width = 8,
             height = 8,
             dpi = 600)

    }




  }


}



#' Flip and/or centre gggenes data
#'
#' Flip and/or centre gggenes data for clearer visualisation
#'
#' @param gggenes_df gggenes data frame (with each gene cluster as "gene cluster")
#' @param gene_on Gene to centre on / use for direction
#' @param direction Direction for gene_on to be in [1/-1] [default = 1] (1 for right, -1 for left)
#' @param length Either "use_max" (use the max gene co-ordinate), or "fasta_length" -
#' @param centre Centre genes relative to all central ("gene_on") genes for each gene cluster in the dataframe [T/F]
#'
#' @return Data frame with each gene cluster co-ordinates set on the centre_gene and re-orientated for uniformity
#' @examples
#' centre_gggenes(gggenes_df = dataframe, gene_on = "TssD", direction = -1, centre = T)


flip_gggenes <- function(gggenes_df, gene_on, direction = 1, length = "use_max", centre = F) {


  # check input

  if (! length %in% c("use_max","fasta_length")) {

    stop("must use either 'use_max' or 'fasta_length")

  }

  if (length == "fasta_length" & ! "cluster_length" %in% colnames(gggenes_df)) {
    stop("lengths of the gene cluster must be labelled as `cluster_length` in the dataframe")
  }

  if (! gene_on %in% gggenes_df$gene) {
    stop("Central gene must be in the `gene` column in the input data frame")
  }



  average_position <- as.integer(mean(subset.data.frame(gggenes_df, gene == gene_on)$start)) # get average starting position for centre gene


  df_list <- list()
  counter <- 0

  for (cluster in unique(gggenes_df$gene_cluster)) {

    counter <- counter + 1
    temp_df <- subset.data.frame(gggenes_df, gene_cluster == cluster)
    gene_direction <- subset.data.frame(temp_df, gene == gene_on)$direction
    if (length(gene_direction) > 1) {
      print(glue::glue("More than one {gene_on}"))
    }
    #print(gene_direction)
    #
    if (nrow(subset.data.frame(temp_df, gene == gene_on)) > 0 ) {
      if (gene_direction == (direction * -1) )  {
        ##### need to reverse the numbers and the directions
        if (length == "use_max") {
          max_length <- max(temp_df$end)
        } else if (length == "fasta_length") {
          max_length <- max(temp_df$cluster_length)
        }

        ## minus the length from the start and stop
        temp_df$start <- max_length - temp_df$start
        temp_df$end <- max_length - temp_df$end
        ## change colnames to switch the start and the stop over
        temp_df <- temp_df %>%
          rename(start = end, end = start) %>%
          mutate(flipped = TRUE)
        #colnames(temp_df) <- c("operon","number","end","start","gene","strand","direction")
        ## change direction:
        temp_df$direction <- temp_df$direction * -1
      } else { # if the gene cluster was already in the correct orientation
        temp_df <- temp_df %>%
          mutate(flipped = FALSE)
      }

      ###### correct the position of each of these if want to centre:
      if (isTRUE(centre)) {
        difference <- average_position + subset.data.frame(temp_df, gene == gene_on)$start
        temp_df$start <- temp_df$start - difference
        temp_df$end <- temp_df$end - difference
      }

    } else { # if couldn't find the gene to centre the "flipping" on
      temp_df <- temp_df %>%
        mutate(flipped = FALSE)
    }
    ## concatenate to list
    df_list[[counter]] <- temp_df
  }


  return(dplyr::bind_rows(df_list))
}








