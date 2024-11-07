

#' Read in annotation from gff file into a tidy dataframe
#'
#' Read in annotation from gff file into a tidy dataframe (including all fields in the decription [9th] column in gff files)
#'
#' @param gff_file Gff file
#' @return dataframe with annotation from gff file
#' @examples
#' read_gff("input.gff")



read_gff <- function(gff_file) {

  gff.data <- read.delim(file = gff_file,
                         sep = "\t",
                         header = F,
                         fill = T,
                         comment.char = "",
                         col.names = c("contig","method","feature","start","end","score","strand","frame","desc")) %>%
    filter(desc != "") %>%
    #extend the description column into multiple columns
    mutate(desc = strsplit(as.character(desc), ";")) %>%
    tidyr::unnest(desc) %>%
    #separate by the "=" character to get new col names, then pivot wider
    tidyr::separate(col = desc, sep = "=", into = c("desc","value"), extra = "merge") %>%
    tidyr::pivot_wider(names_from = "desc", values_from = "value", values_fn = list) %>%
    #mutate across all so that lists are removed and they're just characters
    mutate(across(everything(), as.character)) %>%
    #but change the start and end back to integers
    mutate(start = as.integer(start),
           end = as.integer(end))


  return(gff.data)

}




#' Reverse a gff file
#'
#' Reverse a gff file and write out to file
#'
#' @param input_gff Gff file
#' @param output_gff Gff file
#' @examples
#' reverse_gff("path/to/input.gff", "path/to/output.gff")

reverse_gff <- function(input_gff, output_gff) {
  # Read the GFF file
  gff_data <- readLines(input_gff)

  # Initialize variables
  features <- list()
  sequence_lines <- c()
  headers <- c()
  fasta <- FALSE
  length <- 0
  fasta_seq <- c()

  for (line in gff_data) {
    if (fasta) {
      sequence_lines <- c(sequence_lines, line)
      next
    }
    if (startsWith(line, "##FASTA")) {
      fasta <- TRUE
      next
    }
    if (startsWith(line, "##sequence-region")) {
      headers <- c(headers, line)
      length <- as.integer(strsplit(line, " ")[[1]][4])
      seq.name <- as.character(strsplit(line, " ")[[1]][2])
      next
    }
    if (startsWith(line, "#")) {
      headers <- c(headers, line)
      next
    }

    toks <- strsplit(line, "\t")[[1]]
    if (!grepl("ID=", toks[9])) {
      next
    }

    name <- strsplit(strsplit(toks[9], "ID=")[[1]][2], ";")[[1]][1]
    features[[name]] <- list(
      contig = toks[1],
      method = toks[2],
      feature = toks[3],
      start = as.integer(toks[4]),
      stop = as.integer(toks[5]),
      score = toks[6],
      strand = ifelse(toks[7] == "+", "-", ifelse(toks[7] == "-", "+", ".")),
      frame = toks[8],
      desc = toks[9]
    )
    features[[name]]$temp_start <- length - features[[name]]$stop + 2
    features[[name]]$stop <- length - features[[name]]$start + 2
    features[[name]]$start <- features[[name]]$temp_start

  }



  # Reverse complement the sequence
  joined_sequence <- paste(sequence_lines[-1], collapse = "")
  reverse_complement_seq <- Biostrings::DNAStringSet(
    Biostrings::reverseComplement(Biostrings::DNAString(joined_sequence))
  )
  #set a name to this?
  names(reverse_complement_seq) <- seq.name
  # reverse_complement_seq <- chartr("ACGTacgt", "TGCAtgca", sapply(strsplit(joined_sequence, NULL)[[1]], rev))
  # reverse_complement_seq <- paste(reverse_complement_seq, collapse = "")

  new_sequence <- c(headers)

  for (name in names(features)) {
    details <- features[[name]]
    new_sequence <- c(new_sequence, sprintf(
      "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s",
      details$contig, details$method, details$feature, details$start,
      details$stop, details$score, details$strand, details$frame, details$desc
    ))
  }

  new_sequence <- c(new_sequence, "##FASTA")


  # Write the new GFF file
  writeLines(new_sequence, output_gff)

  #then write the fasta
  Biostrings::writeXStringSet(x = reverse_complement_seq,
                              filepath = output_gff,
                              append = T)


}


